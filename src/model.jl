"""
    StageProblem(scenarios)

Constructs a subproblem with the StageDataExt extension
"""
function StageProblem(scenarios::Int=1)
    sp = Model()
    sp.ext[:data] = StageDataExt(scenarios)
    return sp
end

"""
    stagedata(subproblem)

Get the StageDataExt from a subproblem
"""
function stagedata(m::Model)
    @assert haskey(m.ext, :data)
    return m.ext[:data]::StageDataExt
end

"""
    issubproblem(model)

Checks if a JuMP model is a SDDP subproblem
"""
issubproblem(m::JuMP.Model) = isa(stagedata(m), StageDataExt)
issubproblem(m) = false

# Some helper functions
JuMP.getlowerbound(c::ConstraintRef) = c.m.linconstr[c.idx].lb
JuMP.getupperbound(c::ConstraintRef) = c.m.linconstr[c.idx].ub

# get the number of forward passes conducted
getn(m::SDDPModel) = m.forwardstorage.n
getmarkov(m::SDDPModel, pass, t) = m.forwardstorage.W[pass][t]
getx(m::SDDPModel, pass, t, idx) = m.forwardstorage.x[pass][idx,t]
getx(m::SDDPModel, pass, t) = m.forwardstorage.x[pass][:,t]
getobj(m::SDDPModel) = m.forwardstorage.obj

function setn!(m::SDDPModel, n::Int)
    m.forwardstorage.n = n
end

# Store the markov state in the forward pass datastructure
function savemarkov!(m::SDDPModel, pass, stage, markov)
    m.forwardstorage.W[pass][stage] = markov
end

# save the objective value in the forward pass datastructure
function saveobj!(m::SDDPModel, pass, obj)
    m.forwardstorage.obj[pass] += obj
end

# save the state solution in the forward pass datastructure
function savex!(m::SDDPModel, sp::Model, pass, stage)
    for (i, var) in enumerate(stagedata(sp).state_vars)
        m.forwardstorage.x[pass][i, stage] = getvalue(var)
    end
end

# TODO: This isn't a proper copy
# function Base.copy(sd::StageDataExt)
#     StageDataExt(
#         sd.state_vars,
#         sd.dual_constraints,
#         sd.theta,
#         sd.last_markov, sd.scenario_constraints, sd.scenario_constraint_names, sd.objective_value, sd.dual_values, sd.stage_profit, sd.value_to_go_bound, sd.beta_quantile, sd.lambda_weight, sd.weightings_matrix, sd.regularisecoefficient, sd.regularisecons, sd.regularisepen)
# end
function Base.copy{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM})
    SDDPModel{T, M, S, X, TM}(
        deepcopy(m.stage_problems),
        copy(m.transition),
        m.initial_markov_state,
        m.scenario_probability,
        m.confidence_interval,
        m.valid_bound,
        m.build_function!,
        deepcopy(m.stagecuts),
        m.solver,
        m.valuetogobound,
        m.forwardstorage,
        m.roundingaccuracy
    )
end

"""
    getstagevalue(subproblem)

This function gets the stage costs accrued in the current stage (i.e. the objective less the cost-to-go).
"""
function getstagevalue(sp::Model)
    @assert issubproblem(sp)
    if stagedata(sp).theta != nothing
        return (getobjectivevalue(sp) - getvalue(stagedata(sp).theta))::Float64
    else
        return getobjectivevalue(sp)::Float64
    end
end

"""
    m = SDDPModel(kwargs...) do subproblem, stage, markovstate

    end

Keyword Arguments:
  initial_markov_state      the scenario at time 0. Model transitions at the start of time period 1
  markov_states             the number of markov states
  scenarios                 the number of stagewise independent scenarios in a markov state
  scenario_probability      vector with probability of each scenario occuring
  sense                     :Max or :Min
  solver                    AbstractMathProgBase solver capable of returning dual variables. Defaults   to Clp.
  stages                    the number of stages
  transition                the transition matrix. Can be N*N where N is the number of scenarios, or M*N*N where M is the number of stages.
  value_to_go_bound         Initial bound on value to go
  rounding_accuracy         Number of decimal places to round cut coefficients to. Use if you run into scaling errors
"""
SDDPModel(
    build_subproblem!::Function;
    initial_markov_state = 0,
    markov_states        = 1,
    scenarios            = 1,
    scenario_probability = nothing,
    sense                = :Max,
    stages               = 1,
    transition           = nothing,
    solver               = ClpSolver(),
    value_to_go_bound    = NaN,
    rounding_accuracy    = 6
    ) = SDDPModel(build_subproblem!, initial_markov_state, markov_states,
            scenarios, scenario_probability, sense, stages, transition,
            solver, value_to_go_bound, rounding_accuracy
        )

function SDDPModel(
        build_subproblem!::Function,
        initial_markov_state::Int,
        markov_states::Int,
        scenarios::Int,
        scenario_probability,
        sense::Symbol,
        stages::Int,
        transition,
        solver::MathProgBase.AbstractMathProgSolver,
        value_to_go_bound,
        rounding_accuracy::Int
        )
    # Some sanity checks
    isnan(value_to_go_bound) && error("You must specify the option [value_to_go_bound] when creating an SDDPModel.")
    @assert stages >= 1
    @assert markov_states >= 1
    @assert scenarios >= 1
    !(sense in [:Max, :Min]) && error("Sense must be either :Max or :Min. Currently = $(sense).")

    transition_matrix, T = maketransition(transition, stages, markov_states)

    if scenario_probability == nothing
        scenario_probability = ones(scenarios) / scenarios
    else
        @assert abs(sum(scenario_probability) - 1) < 1e-5
    end

    my_inf = (sense==:Max?Inf:-Inf)
    sense_type = (sense==:Max?Max:Min)

    m = SDDPModel{stages, markov_states, scenarios,sense_type, T}(
        Array(JuMP.Model, (stages, markov_states)),
        transition_matrix,
        initial_markov_state,
        WeightVec(scenario_probability),
        (-my_inf, -my_inf),
        my_inf,
        build_subproblem!,
        Array(StageCuts, (stages, markov_states)),
        solver,
        value_to_go_bound,
        ForwardPassData(),
        rounding_accuracy
        )

    create_subproblems!(m)
    initialise_cutstorage!(m)

    return m
end

function initialise_cutstorage!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM})
    for t=1:T
        for i=1:M
            m.stagecuts[t,i] = StageCuts(subproblem(m, t, i), m.valuetogobound)
        end
    end
end

function maketransition(transition, stages, markovstates)
    if transition==nothing
        # User hasn't specified transition matrix. Assume uniform
        return 1 / markovstates * ones(markovstates, markovstates), 2
    end
    T = length(size(transition))
    if T==2
        # Stage independent transition probabilities
        @assert is_valid_transition_matrix(transition, markovstates)
        return transition, 2
    elseif T==1
        # Vector of transition matrices
        @assert length(transition) == stages
        transitionout = Array(Float64, (markovstates, markovstates, stages))
        for t=1:stages
            @assert is_valid_transition_matrix(transition[t], markovstates)
            for i=1:markovstates
                for j=1:markovstates
                    transitionout[i, j, t] = transition[t][i,j]
                end
            end
        end
        return transitionout, 3
    end
    error("Transition matrix was wrong shape.")
end

"""
This fuction checks that the matrix T is a valid markov transition matrix
"""
function is_valid_transition_matrix(T::Array{Float64, 2}, n::Int)
    if size(T)[1] != size(T)[2]
        error("Transition matrix must be square")
    end
    @assert size(T)[1] == n
    for i=1:n
        if abs(sum(T[i,:]) - 1.) > 1e-8
            return false
        end
    end
    return true
end


# This will probably break at some point
# https://groups.google.com/forum/#!searchin/julia-users/arguments$20anonymous$20function/julia-users/QcgdNZd-sI8/dgIYSxozZzUJ
if VERSION > v"0.5-"
    arglength(f::Function)=length(Base.uncompressed_ast(methods(f).defs.func).args[1])-1
else
    arglength(f::Function)=length(Base.uncompressed_ast(f.code.def).args[1])
end

"""
    create_subproblems!(SDDPModel, solver, valuebound)

Create all the JuMP subproblems for the SDDPModel
"""
function create_subproblems!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM})
    # For every stage
    for stage=1:T
        # For every markov state
        for markov_state=1:M
            # Create a new stage problem with [scenarios] number of scenarios
            sp = StageProblem(S)
            # Set the solver
            setsolver(sp, m.solver)
            if arglength(m.build_function!)==3
                # User has specified (subproblem, stage, markov states)
                m.build_function!(sp, stage, markov_state)
            elseif arglength(m.build_function!)==2
                # No markov states specified
                # Double check they only mean one markov state
                @assert M == 1
                m.build_function!(sp, stage)
            else
                error("Invalid number of arguments")
            end
            # Initialise storage for scenario duals now we know the number of state variables
            for i in 1:S
                push!(stagedata(sp).dual_values, zeros(length(stagedata(sp).state_vars)))
            end

            # Initialise regularisation variables
            stagedata(sp).regularisepen = @variable(sp, regularisepen >= 0)
            @constraints(sp, begin
                regpos, stagedata(sp).regularisepen + sum{v, v in stagedata(sp).state_vars} >= 0
                regneg, stagedata(sp).regularisepen - sum{v, v in stagedata(sp).state_vars} >= 0
            end)
            push!(stagedata(sp).regularisecons, regpos)
            push!(stagedata(sp).regularisecons, regneg)

            if stage<T
                addtheta!(X, sp, m.valuetogobound)
            end
            set_objective!(X, sp)


            stagedata(sp).weightings_matrix = zeros(M,S)
            for j=1:M
                for s=1:S
                    stagedata(sp).weightings_matrix[j,s] = transitionprobability(m, stage, markov_state, j) * m.scenario_probability[s]
                end
            end

            # Store the stage problem
            m.stage_problems[stage, markov_state] = sp
        end
    end
    return
end

# Add the value/cost to go variable depending on model sense
addtheta!(::Type{Min}, sp::Model, value_to_go_bound) = (stagedata(sp).theta = @variable(sp, theta >= value_to_go_bound))
addtheta!(::Type{Max}, sp::Model, value_to_go_bound) = (stagedata(sp).theta = @variable(sp, theta <= value_to_go_bound))

# Set the objective depending on model sense
setobj!(::Type{Min}, sp, aff) = @objective(sp, Min, aff)
setobj!(::Type{Max}, sp, aff) = @objective(sp, Max, aff)

# Linear regularisation term
function regularise!(regularisation::LinearRegularisation, sense::Sense, sp)
    stagedata(sp).regularisecoefficient *= regularisation.decayrate
    return sense!(sense, stagedata(sp).regularisecoefficient * stagedata(sp).regularisepen)
end

function oldvalue(v)
    ov = getvalue(v)
    if isnan(ov) || ov == -Inf || ov == Inf
        ov = 0.
    end
    ov
end

# Quadratic regularisation term
function regularise!(regularisation::QuadraticRegularisation, sense::Sense, sp)
    stagedata(sp).regularisecoefficient *= regularisation.decayrate
    @expression(sp, regulariser, sum{(v - oldvalue(v))^2, v in stagedata(sp).state_vars})
    return sense!(sense, stagedata(sp).regularisecoefficient * regulariser)
end

# Dummy regularisation
regularise!(ty, sense::Sense, sp) = 0.

# Get regularisation sense correct
sense!(::Type{Min}, expr) = expr
sense!(::Type{Max}, expr) = -expr

# Return future cost
#    TODO: this isn't type stable
function futurecost!(sp)
    if stagedata(sp).theta != nothing
        return stagedata(sp).theta
    end
    return 0.
end

# Set potentially regularised objective
function set_objective!(regularise, sense::Sense, sp::Model)
    aff = stagedata(sp).stage_profit + sense!(sense, stagedata(sp).relaxationpenalties)
    aff += futurecost!(sp)
    aff += regularise!(regularise, sense, sp)
    setobj!(sense, sp, aff)
end
# No regularisation specified
set_objective!(sense::Sense, sp::Model) = set_objective!(nothing, sense, sp)

set_nonregularised_objective!(regularisation::Regularisation, sense::Sense, sp::Model) = set_objective!(nothing, sense, sp)
set_nonregularised_objective!(::NoRegularisation, sense::Sense, sp::Model) = nothing

set_regularised_objective!(regularisation::Regularisation, sense::Sense, sp::Model) = set_objective!(regularisation, sense, sp)
set_regularised_objective!(regularisation::NoRegularisation, sense::Sense, sp::Model) = nothing

# Get markov transition. r is used for antithetic variates
function transition{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t::Int, i::Int, r::Float64)
    for j=1:M
        k = transitionprobability(m, t, i, j)
        if r <= k
            return j
        else
            r -= k
        end
    end
    return M
end
transition{T, S, X, TM}(m::SDDPModel{T, 1, S, X, TM}, t::Int, i::Int, r::Float64) = 1
transition{T, S, X, TM}(m::SDDPModel{T, 1, S, X, TM}, t::Int, i::Int) = 1
transition{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t::Int, i::Int) = transition(m, t, i, rand())
function transitionprobability{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t::Int, i::Int, j::Int)
    if i==0
        return 1./M
    else
        return transitionprobability(m.transition, t, i, j)::Float64
    end
end
transitionprobability(TM::Array{Float64, 2}, t, i, j) = TM[i,j]
transitionprobability(TM::Array{Float64, 3}, t, i, j) = TM[i,j, t]

function setbound!(::Type{Min}, m::SDDPModel, obj)
    if obj > getBound(m)
        m.valid_bound = obj
    end
end
function setbound!(::Type{Max}, m::SDDPModel, obj)
    if obj < getBound(m)
        m.valid_bound = obj
    end
end
setbound!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, obj) = setbound!(X, m, obj)

"""
    setbound!(SDDPModel)

Calculate the upper bound of the first stage problem
"""
function setbound!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM})
    # Initialise
    obj = 0.0
    if m.initial_markov_state == 0
        # Lets average over all first stage probles (markov states x scenarios)
        for i=1:M
            sp = subproblem(m, 1, i)        # get subproblem
            for s=1:S
                load_scenario!(sp, s)   # realise scenario
                forwardsolve!(sp)          # solve
                obj += transitionprobability(m, 1, m.initial_markov_state, i) * m.scenario_probability[s] * getobjectivevalue(sp)
            end
        end
    else
        # Lets just  average over the scenarios in the initial markov state
        sp = subproblem(m, 1, m.initial_markov_state) # get subproblem
        for s=1:S
            load_scenario!(sp, s)                  # load scenario
            forwardsolve!(sp)                         # solve
            obj += m.scenario_probability[s] * getobjectivevalue(sp)
        end
    end
    setbound!(m, obj)    # Update the bound
end

function estimatebound(obj::Vector{Float64}, conflevel)
    if length(obj) > 1
        return confidenceinterval(obj, conflevel)
    else
        return (mean(obj), mean(obj))
    end
end

# Construct a confidence interval
function confidenceinterval(x, conf_level=0.95)
    tstar = quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    SE = std(x)/sqrt(length(x))
    lo, hi = mean(x) + [-1, 1] * tstar * SE
    return (lo, hi)
end
setconfidenceinterval!{T<:Real}(m::SDDPModel, v::Tuple{T,T}) = (m.confidence_interval = v)

function setbound!(log::SolutionLog, m::SDDPModel)
    setbound!(m)
    log.bound = m.valid_bound
    return
end

# a wrapper for estimating the confidence interval
function setconfidenceinterval!(log, m, conflevel)
    setconfidenceinterval!(m, estimatebound(getobj(m), conflevel))
    log.ci_lower = log.ci_upper = mean(m.confidence_interval)
    return
end

isconverged(::Type{Min}, ci::NTuple{2, Float64}, bound::Float64) = ci[1] < bound
isconverged(::Type{Max}, ci::NTuple{2, Float64}, bound::Float64) = ci[2] > bound
isconverged{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, ci, bound) = isconverged(X, ci, bound)

function solve!(sp::Model)
    @assert issubproblem(sp)
    status = solve(sp)
    # Catch case where we aren't optimal
    if status != :Optimal
        warn("SDDP subproblem not optimal (stats=$(status)). Assuming numerical infeasibility so rebuilding model from stored cuts.")
        sp.internalModelLoaded = false
        status = solve(sp)
        if status != :Optimal
            JuMP.writeMPS(sp, "subproblem_proc$(myid())_$(string(hash(now()))).mps")
            error("SDDP Subproblems must be feasible. Current status: $(status). I tried rebuilding from the JuMP model but it didn't work so I wrote you an MPS file.")
        end
    end
end
