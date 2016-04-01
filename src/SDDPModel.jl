type StageData
    state_vars::Vector{JuMP.Variable}
    dual_constraints::Vector{JuMP.ConstraintRef}
    theta::Union{Void, JuMP.Variable}
    old_scenario::Tuple{Int, Int}
    scenario_constraints::Vector{Tuple{Any, Vector{Any}}}
    last_scenario::Int
    current_scenario::Int
    objective_value::Vector{Float64}
    dual_values::Vector{Vector{Float64}}
    stage_profit
end
StageData(scenarios::Int=1) = StageData(Variable[], ConstraintRef[], nothing, (0,0), Tuple{Any, Vector{Any}}[], 0, 0, zeros(scenarios), Vector{Float64}[], nothing)

"""
Instaniates a new StageProblem which is a JuMP.Model object with an extension
dictionary.
"""
function StageProblem(scenarios::Int=1)
    sp = Model()
    sp.ext[:data] = StageData(scenarios)
    return sp
end

function stagedata(m::Model)
    @assert haskey(m.ext, :data)
    return m.ext[:data]::StageData
end

# Check this is a StageProblem
is_sp(m::JuMP.Model) = isa(stagedata(m), StageData)
is_sp(m) = false

type SDDPModel{M,N,S,T}
    sense::Symbol
    stage_problems::Array{JuMP.Model}
    transition::Array{Union{Float64, Array{Float64, 2}}, T}
    init_markov_state::Int
    confidence_interval::Tuple{Float64, Float64}
    valid_bound::Float64
    QUANTILE::Float64
    LPSOLVER::MathProgBase.AbstractMathProgSolver
    value_to_go_bound::Float64
    beta_quantile::Float64
    risk_lambda::Float64
    weightings_matrix::Array{Float64, 2}
    cuts_filename::Union{Void, ASCIIString}
end

"""
Creates an SDDPModel type

    SDDPModel(kwargs...)

Keyword Arguments:
conf_level        - confidence level for the lower bound
cuts_filename     - ASCIIString filename for cut output. If specified cuts will be written.
initial_markov_state - the scenario at time 0. Model transitions at the start of time period 1
markov_states     - the number of markov states
scenarios         - the number of stagewise independent scenarios in a markov state
sense             - :Max or :Min
solver            - AbstractMathProgBase solver capable of returning dual variables. Defaults to Clp.
stages            - the number of stages
transition        - the transition matrix. Can be N*N where N is the number of scenarios, or M*N*N where M is the number of stages.
value_to_go_bound - Initial bound on value to go
"""
function SDDPModel(;
    conf_level=0.95,
    cuts_filename=nothing,
    initial_markov_state=0,
    markov_states=1,
    scenarios=1,
    sense=:Max,
    solver=ClpSolver(),
    stages=1,
    transition=nothing,
    value_to_go_bound=1000
    )

    # Check non-zero stages, markov states and scenarios
    @assert stages >= 1
    @assert markov_states >= 1
    @assert scenarios >= 1

    if !(sense in [:Max, :Min])
        error("Sense must be either :Max or :Min. Currently = $(sense).")
    end

    if transition==nothing
        # User hasn't specified transition matrix. Assume uniform
        transition = 1 / markov_states * ones(markov_states, markov_states)
    end
    T = length(size(transition))
    if T==2
        # Stage independent transition probabilities
        @assert is_valid_transition_matrix(transition, markov_states)
    elseif T==1
        # Vector of transition matrices
        @assert length(transition) == stages
        for i=1:stages
            @assert is_valid_transition_matrix(transition[i], markov_states)
        end
    end
    # TODO :: Case where Transition matrix is Array{Float64, 3}

    my_inf = (sense==:Max?Inf:-Inf)
    SDDPModel{stages,markov_states,scenarios,T}(sense, Array(JuMP.Model, (stages, markov_states)), transition, initial_markov_state, (-my_inf, -my_inf), my_inf, conf_level, solver, value_to_go_bound, 1., 1., zeros(markov_states, scenarios), cuts_filename)
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
        if abs(sum(T[i,:]) - 1.) > 1e-10
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

function is_zero_objective{T}(ex::JuMP.GenericQuadExpr{T, JuMP.Variable})
    return length(ex.qcoeffs) == 0 && length(ex.qvars1) == 0 && length(ex.qvars2) == 0 && is_zero_objective(ex.aff)
end
function is_zero_objective{T}(ex::JuMP.GenericAffExpr{T, JuMP.Variable})
    return length(ex.coeffs) == 0 && length(ex.vars) == 0 && ex.constant == 0.0
end

function SDDPModel(
    build_subproblem!::Function;
    conf_level=0.95,
    cuts_filename=nothing,
    initial_markov_state=0,
    markov_states=1,
    scenarios=1,
    sense=:Max,
    stages=1,
    transition=nothing,
    solver=ClpSolver(),
    value_to_go_bound=NaN
    )
    if isnan(value_to_go_bound)
        error("You must specify the option [value_to_go_bound] when creating an SDDPModel.")
    end

    # Initialise SDDPModel object
    m = SDDPModel(
        sense=sense,
        stages=stages,
        markov_states=markov_states,
        scenarios=scenarios,
        transition=transition,
        initial_markov_state=initial_markov_state,
        conf_level=conf_level,
        solver=solver,
        value_to_go_bound=value_to_go_bound,
        cuts_filename=cuts_filename
    )

    # For every stage
    for stage=1:stages
        # For every markov state
        for markov_state=1:markov_states
            # Create a new stage problem with [scenarios] number of scenarios
            sp = StageProblem(scenarios)

            # Set the solver
            setSolver(sp, m.LPSOLVER)

            # User has specified markov states
            if arglength(build_subproblem!)==3
                build_subproblem!(sp, stage, markov_state)

            # No markov states specified
            elseif arglength(build_subproblem!)==2
                # Double check they only mean one markov state
                @assert markov_states == 1

                build_subproblem!(sp, stage)
            else
                error("Invalid number of arguments")
            end

            # Initialise storage for scenario duals now we know the number of state variables
            for i in 1:length(stagedata(sp).state_vars)
                push!(stagedata(sp).dual_values, zeros(scenarios))
            end

            # If the user hasn't specified an objective
            if is_zero_objective(getObjective(sp))
                if stage==stages
                    # If its the last stage then its just the stage profit
                    @setObjective(sp, sense, stagedata(sp).stage_profit)
                else
                    # Otherwise create a value/cost to go variable
                    if sense==:Max
                        @defValueToGo(sp, theta <= value_to_go_bound)
                    else
                        @defValueToGo(sp, theta >= value_to_go_bound)
                    end
                    # Set the objective
                    @setObjective(sp, sense, stagedata(sp).stage_profit + theta)
                end
            end
            # Store the stage problem
            m.stage_problems[stage, markov_state] = sp
        end
    end

    return m
end

"""
So we can copy an SDDPModel
"""
function Base.copy{M,N,S,T}(m::SDDPModel{M,N,S,T})
    SDDPModel{M,N,S,T}(m.sense, deepcopy(m.stage_problems), copy(m.transition), m.init_markov_state, m.confidence_interval, m.valid_bound, m.QUANTILE, m.LPSOLVER, m.value_to_go_bound, m.beta_quantile, m.risk_lambda, m.weightings_matrix, m.cuts_filename)
end

"""
Transition markov states

Inputs:
m            - The SDDPModel object
stage        - the current stage
markov_state - the current markov state

Returns:
    the new markov state::Int
"""
function transition{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, markov_state::Int)
    r = rand(Float64)
    for new_markov_state=1:N
        k = get_transition(m, stage, markov_state, new_markov_state)
        if r <= k
            return new_markov_state
        else
            r -= k
        end
    end
    return N
end
function get_transition{M,N,S}(m::SDDPModel{M,N,S,2}, stage::Int, current_markov_state::Int, new_markov_state::Int)
    if current_markov_state==0
        return 1./N
    else
        return m.transition[current_markov_state,new_markov_state]::Float64
    end
end
function get_transition{M,N,S}(m::SDDPModel{M,N,S,1}, stage::Int, current_markov_state::Int, new_markov_state::Int)
    if current_markov_state==0
        1./N
    else
        m.transition[stage][current_markov_state,new_markov_state]::Float64
    end
end

"""
This function loads cuts from a file.
"""
function load_cuts!{M,N,S,T}(m::SDDPModel{M,N,S,T}, filename::ASCIIString)
    open(filename, "r") do f
        while true
            line = readline(f)
            (line == nothing || line == "") && break
            line = split(strip(line), ",")
            @assert length(line) >= 3
            stage = parse(Int, line[1])
            markov_state = parse(Int, line[2])
            sp = m.stage_problems[stage, markov_state]
            theta = parse(Float64, line[3])

            @assert length(line) == (3 + length(stagedata(sp).state_vars))
            if length(line) > 4
                xcoeff = map(x->parse(Float64, x), line[4:end])
            else
                xcoeff = [parse(Float64, line[4])]
            end
            @addConstraint(sp, stagedata(sp).theta <= theta + sum{xcoeff[i]*getVar(sp, v), (i, v) in enumerate(stagedata(sp).state_vars)})
        end
    end
end
function load_cuts!{M,N,S,T}(m::SDDPModel{M,N,S,T})
    if m.cuts_filename == nothing
        error("Please specify a file to load cuts from.")
    end
    load_cuts!(m, m.cuts_filename)
end
