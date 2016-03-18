# TODO at the moment we assume uniform scenario probability in each markov state

module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp# , Gurobi
using Formatting
using Distributions

export SDDPModel,
    @defStateVar, @defValueToGo, @addScenarioConstraint, @setStageProfit,
    simulate # addStageProblem!,

type SDDPModel{M,N,S,T}
    sense::Symbol
    stage_problems::Array{JuMP.Model}
    transition::Array{Union{Float64, Array{Float64, 2}}, T}
    init_markov_state::Int
    confidence_interval::Tuple{Float64, Float64}
    valid_bound::Float64
    status::Symbol
    QUANTILE::Float64
    LPSOLVER::MathProgBase.AbstractMathProgSolver
    theta_bound::Float64
    beta_quantile::Float64
    risk_lambda::Float64
end

include("macros.jl")

"""
Creates an SDDPModel type

    SDDPModel([;sense, stages, markov_states, transition, init_markov_state, conf_level, solver, abs_tol, rel_tol])

Inputs:
sense      - :Max or :Min
stages     - the number of stages
markov_states - the number of scenarios
transition - the transition matrix. Can be N*N where N is the number of scenarios, or M*N*N where M is the number of stages.
init_markov_state - the scenario at time 0. Model transitions at the start of time period 1
conf_level - confidence level for the lower bound
solver   - AbstractMathProgBase solver capable of returning dual variables. Defaults to Clp.
"""
function SDDPModel(;
    sense=:Max,
    stages=1,
    markov_states=1,
    scenarios=1,
    transition=nothing,
    initial_markov_state=0,
    conf_level=0.95,
    solver=ClpSolver(),
    theta_bound=1000)

    @assert stages >= 1
    @assert markov_states >= 1
    @assert scenarios >= 1

    if !(sense in [:Max, :Min])
        error("Sense must be either :Max or :Min. Currently = $(sense).")
    end

    if transition==nothing
        # Uniform transition matrix
        transition = 1 / markov_states * ones(markov_states, markov_states)
    end
    T = length(size(transition))
    # Check valid transition matrix

    if T==2
        @assert size(transition)[1] == size(transition)[2] == markov_states
        for i=1:markov_states
            @assert abs(sum(transition[i,:]) - 1.) < 1e-10
        end
    elseif T==1
        @assert length(transition) == stages
        for i=1:stages
            @assert size(transition[i])[1] == size(transition[i])[2] == markov_states
            for j=1:markov_states
                @assert abs(sum(transition[i][j,:]) - 1.) < 1e-10
            end
        end
    end
    my_inf = (sense==:Max?Inf:-Inf)
    SDDPModel{stages,markov_states,scenarios,T}(sense, Array(JuMP.Model, (stages, markov_states)), transition, initial_markov_state, (-my_inf, -my_inf), my_inf, :Unconverged, conf_level, solver, theta_bound, 1., 1.)
end

# This will probably break at some point
# https://groups.google.com/forum/#!searchin/julia-users/arguments$20anonymous$20function/julia-users/QcgdNZd-sI8/dgIYSxozZzUJ
if VERSION > v"0.5-"
    arglength(f::Function)=length(Base.uncompressed_ast(methods(f).defs.func).args[1])-1
else
    arglength(f::Function)=length(Base.uncompressed_ast(f.code.def).args[1])
end

function is_zero_objective{T}(ex::JuMP.GenericQuadExpr{T, JuMP.Variable})
    return length(ex.qcoeffs) == 0 &&
        length(ex.qvars1) == 0 &&
        length(ex.qvars2) == 0 &&
        is_zero_objective(ex.aff)
end
function is_zero_objective{T}(ex::JuMP.GenericAffExpr{T, JuMP.Variable})
    return length(ex.coeffs) == 0 &&
    length(ex.vars) == 0 &&
    ex.constant == 0.0
end

function SDDPModel(build_subproblem!::Function;
    sense=:Max,
    stages=1,
    markov_states=1,
    scenarios=1,
    transition=nothing,
    initial_markov_state=0,
    conf_level=0.95,
    solver=ClpSolver(),
    theta_bound=1e3)

    m = SDDPModel(sense=sense,
    stages=stages,
    markov_states=markov_states,
    scenarios=scenarios,
    transition=transition,
    initial_markov_state=initial_markov_state,
    conf_level=conf_level,
    solver=solver,
    theta_bound=theta_bound)

    for stage=1:stages
        for markov_state=1:markov_states
            # for scenario=1:scenarios
            sp = StageProblem(scenarios)
            setSolver(sp, m.LPSOLVER)
            if arglength(build_subproblem!)==3
                build_subproblem!(sp, stage, markov_state)
            elseif arglength(build_subproblem!)==2
                @assert markov_states == 1
                build_subproblem!(sp, stage)
            else
                error("Invalid number of arguments")
            end

            for v in sp.ext[:state_vars]
                sp.ext[:DualValues][v] = zeros(scenarios)
            end

            if is_zero_objective(getObjective(sp))
                if stage==stages
                    @setObjective(sp, sense, sp.ext[:StageProfit])
                else
                    if sense==:Max
                        @defValueToGo(sp, theta <= theta_bound)
                    else
                        @defValueToGo(sp, theta >= theta_bound)
                    end
                    @setObjective(sp, sense, sp.ext[:StageProfit] + theta)
                end
            end
            m.stage_problems[stage, markov_state] = sp
            # end
        end
    end

    return m
end

function Base.copy{M,N,S,T}(m::SDDPModel{M,N,S,T})
    SDDPModel{M,N,S,T}(m.sense, deepcopy(m.stage_problems), copy(m.transition), m.init_markov_state, m.confidence_interval, m.valid_bound, m.status, m.QUANTILE, m.LPSOLVER, m.theta_bound, m.beta_quantile, m.risk_lambda)
end

function transition{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, scenario::Int)
    r = rand(Float64)
    for i=1:N
        k = get_transition(m, stage, scenario, i)
        if r <= k
            return i
        else
            r -= k
        end
    end
end
function get_transition{M,N,S}(m::SDDPModel{M,N,S,2}, stage::Int, s1::Int, s2::Int)
    if s1==0
        return 1/N
    else
        return m.transition[s1,s2]
    end
end
function get_transition{M,N,S}(m::SDDPModel{M,N,S,1}, stage::Int, s1::Int, s2::Int)
    s1==0?1/N : m.transition[stage][s1,s2]
end

"""
Instaniates a new StageProblem which is a JuMP.Model object with an extension
dictionary.
"""
function StageProblem(scenarios::Int=1)
    sp = Model()
    sp.ext[:is_sp] = true
    sp.ext[:state_vars] = Symbol[]
    sp.ext[:duals] = Dict{Symbol, ConstraintRef}()
    sp.ext[:theta] = nothing
    sp.ext[:old_scenario] = (0,0)
    sp.ext[:Scenarios] = Tuple{Any, Vector{Any}}[]
    sp.ext[:LastScenario] = 0
    sp.ext[:CurrentScenario] = 0
    sp.ext[:LastObjectives] = zeros(scenarios)
    sp.ext[:DualValues] = Dict{Symbol, Vector{Float64}}()
    return sp
end

# """
# Add a new stage problem to the SDDPModel.
#
# Usage
#
#     addStageProblem!(m, stage, markov_state) do sp
#         @defVar(sp, x)
#     end
#
# """
# function addStageProblem!{M,N,S,T}(build_model!::Function, m::SDDPModel{M,N,S,T}, stage_idx::Int, scenario_idx::Int)
#     sp = StageProblem()
#     setSolver(sp, m.LPSOLVER)
#     build_model!(sp)
#
#     # remove last theta
#     if stage_idx == M && sp.ext[:theta] != nothing
#         @addConstraint(sp, sp.ext[:theta] == 0)
#     end
#
#     m.stage_problems[stage_idx, scenario_idx] = sp
# end
# addStageProblem!(build_model!::Function, m::SDDPModel, stage_idx::Int) = addStageProblem!(build_model!, m, stage_idx, 1)
#
# Check this is a StageProblem
is_sp(m::JuMP.Model) = haskey(m.ext, :is_sp) && m.ext[:is_sp]
is_sp(m) = false

# """
# A helper function for the JuMP.getValue so it can extract the values from and SDDP problem.
# """
# function JuMP.getValue(m::SDDPModel, v::Symbol)
#     x = Array(Any, length(m.stage_problems))
#     for (idx, sp) in enumerate(m.stage_problems)
#         try
#             x[idx] = getValue(getVar(sp, v))
#         catch
#         end
#     end
#     if issubtype(typeof(x[1]), Number)
#         return convert(Vector{Float64}, x)
#     else
#         # Deal with more complicated values
#         return x
#     end
# end

# Some helper functions
JuMP.getLower(c::ConstraintRef) = c.m.linconstr[c.idx].lb
JuMP.getUpper(c::ConstraintRef) = c.m.linconstr[c.idx].ub


function getCloseCIBound(m::SDDPModel)
    if m.sense==:Max
        m.confidence_interval[2]
    else
        m.confidence_interval[1]
    end
end
function getFarCIBound(m::SDDPModel)
    if m.sense==:Max
        m.confidence_interval[1]
    else
        m.confidence_interval[2]
    end
end
getBound(m::SDDPModel) = m.valid_bound

setCI!{T<:Real}(m::SDDPModel, v::Tuple{T,T}) = (m.confidence_interval = v)
setBound!{T<:Real}(m::SDDPModel, v::T) = (m.valid_bound = v)

function atol(m)
    if m.sense==:Max
        getBound(m) - getCloseCIBound(m)
    else
        getCloseCIBound(m) - getBound(m)
    end
end

function rtol(m)
    abs(getBound(m)) != Inf?atol(m) / abs(getBound(m)):Inf
end

"""
Get the value of the current stage (i.e. not including the value to go)
"""
function get_true_value(m::Model)
    @assert is_sp(m)
    if m.ext[:theta] != nothing
        # TODO: is this correct for minimisation
        return getObjectiveValue(m) - getValue(m.ext[:theta])
    else
        return getObjectiveValue(m)
    end
end

"""
Calculate the upper bound of the first stage problem
"""
function set_valid_bound!{M,N,S,T}(m::SDDPModel{M,N,S,T})
    obj = 0.0
    if m.init_markov_state == 0
        for i=1:N
            for s=1:S
                obj += get_transition(m, 1, m.init_markov_state, i) * m.stage_problems[1, i].ext[:LastObjectives][s]
            end
        end
        obj /= S
    else
        for s=1:S
            obj += m.stage_problems[1, m.init_markov_state].ext[:LastObjectives][s]
        end
        obj /= S
    end
    if m.sense==:Max && obj < getBound(m)
        setBound!(m, obj)
    elseif m.sense==:Min && obj > getBound(m)
        setBound!(m, obj)
    end
end
# function set_valid_bound!{M,N,S,T}(m::SDDPModel{M,N,S,T})
#     obj = 0.0
#     if m.init_markov_state == 0
#
#         P = zeros(S*N)
#         i=1
#         for mkv=1:N
#             for s=1:S
#                 P[i] = get_transition(m, 1, m.init_markov_state, mkv)/S
#                 i+=1
#             end
#         end
#         w = risk_averse_weightings(vcat([sp.ext[:LastObjectives] for sp in m.stage_problems[1, :]]...), P,  m.beta_quantile, m.sense==:Max)
#         i=1
#         for mkv=1:N
#             for s=1:S
#                 obj += w[i] * m.stage_problems[1, mkv].ext[:LastObjectives][s]
#                 i+=1
#             end
#         end
#
#     else
#
#         P = ones(S) / S
#         w = risk_averse_weightings(m.stage_problems[1, m.init_markov_state].ext[:LastObjectives], P,  m.beta_quantile, m.sense==:Max)
#         for s=1:S
#             obj += w[s] * m.stage_problems[1, m.init_markov_state].ext[:LastObjectives][s]
#         end
#
#     end
#
#
#     if m.sense==:Max && obj < getBound(m)
#         setBound!(m, obj)
#     elseif m.sense==:Min && obj > getBound(m)
#         setBound!(m, obj)
#     end
# end

"""
Solve a StageProblem.

    solve!(sp::Model, m::SDDPModel)
"""
function solve!(sp::Model, m::SDDPModel)
    # Catch case where we aren't optimal
    status = solve(sp)
    if status != :Optimal
        # if status == :Infeasible
        #     info("Printing IIS")
        #     grb_model = MathProgBase.getrawsolver(getInternalModel(sp))
        #     Gurobi.computeIIS(grb_model)
        #     num_constrs = Gurobi.num_constrs(grb_model)
        #     iis_constrs = Gurobi.get_intattrarray(grb_model, "IISConstr",  1, num_constrs)
        #     @show sp.linconstr[find(iis_constrs)]
        # elseif status == :Unbounded
        #     @show sp.colVal
        #     @show Variable(sp, length(sp.colVal))
        # end
        error("SDDP Subproblems must be feasible. Current status: $(status).")
    end
    s = sp.ext[:CurrentScenario]
    sp.ext[:LastObjectives][s] = getObjectiveValue(sp)
    for v in sp.ext[:state_vars]
        sp.ext[:DualValues][v][s] = getDual(sp.ext[:duals][v])
    end
    return
end

function load_scenario!(m::Model, scenario::Int)
    m.ext[:LastScenario] = m.ext[:CurrentScenario]
    for (c, Ω) in m.ext[:Scenarios]
        if m.ext[:LastScenario] == 0
            old_scenario = 0.
        else
            old_scenario = Ω[m.ext[:LastScenario]]
        end
        chgConstrRHS(c, getRHS(c) - old_scenario + Ω[scenario])
    end
    m.ext[:CurrentScenario] = scenario
    return
end


"""
Simulate SDDP model and return variable solutions contained in [vars]
"""
function simulate{M,N,S,T}(m::SDDPModel{M,N,S,T}, n::Int, vars::Vector{Symbol})
    results = Dict{Symbol, Any}(:Objective=>zeros(Float64,n))
    for v in vars
        results[v] = Array(Vector{Any}, M)
        for i=1:M
            results[v][i] = Array(Any, n)
        end
    end

    for pass=1:n
        if m.init_markov_state==0
            markov = transition(m, 1, m.init_markov_state)
        else
            markov = m.init_markov_state
        end
        for stage=1:M
            old_markov = markov

            scenario = rand(1:S)

            sp = m.stage_problems[stage, markov]  # corresponding stage problem
            load_scenario!(sp, scenario)

            solve!(sp, m)                           # solve

            results[:Objective][pass] += get_true_value(sp)         # Add objective

            # pass forward if necessary
            if stage < M
                pass_states_forward!(m, stage, old_markov)
                # transition to new scenario
                markov = transition(m, stage, old_markov)
            end

            for v in vars
                results[v][stage][pass] = getValue(getVar(sp, v))
            end

        end
    end
    return results
end

"""
Forward pass of the SDDP algorithm
"""
function forward_pass!{M,N,S,T}(m::SDDPModel{M,N,S,T}, npasses::Int=1, print_info::Bool=false)
    scenario = 0                                    # Initialise variable
    obj = zeros(npasses)                            # Initialise the objective values

    for pass=1:npasses
        if m.init_markov_state==0
            markov = transition(m, 1, m.init_markov_state)
        else
            markov = m.init_markov_state
        end
        for stage=1:M
            old_markov = markov

            scenario = rand(1:S)

            sp = m.stage_problems[stage, markov]  # corresponding stage problem
            load_scenario!(sp, scenario)

            solve!(sp, m)                                   # solve

            obj[pass] += get_true_value(sp)                 # Add objective

            # pass forward if necessary
            if stage < M
                pass_states_forward!(m, stage, old_markov)

                # transition to new scenario
                markov = transition(m, stage, old_markov)
            end
        end
    end

    # set new lower bound
    if abs(m.risk_lambda - 1) < 1e-5
        if npasses > 1
            _obj = t_test(obj, conf_level=m.QUANTILE)
        else
            _obj = (obj[1], obj[1])
        end
    else
        _obj = cvar_estimate(obj, m.beta_quantile, m.risk_lambda)
    end
    if (m.sense==:Max && _obj[1] > getFarCIBound(m)) || (m.sense==:Min && _obj[2] < getFarCIBound(m))
        setCI!(m, _obj)
    end

end

function cvar_estimate(x, beta=1., lambda=1., n=1000)
    y = zeros(100)#div(length(x), 100))
    z = zeros(n)
    for i=1:n
        rand!(y, x)
        z[i] = cvar(y, beta, lambda)
    end
    t_test(z)
end

function cvar{T}(x::Vector{T}, beta::Float64=1., lambda::Float64=1.)
    @assert beta >= 0 && beta <= 1.
    @assert lambda >= 0 && lambda <= 1.
    lambda * mean(x) + (1 - lambda) * mean(x[x.<quantile(x, beta)])
end



function t_test(x::Vector; conf_level=0.95)
    tstar = quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    SE = std(x)/sqrt(length(x))
    lo, hi = mean(x) + [-1, 1] * tstar * SE
    return (lo, hi)#, mean(x))
end

function pass_states_forward!{M,N,S,T}(m::SDDPModel{M,N,S,T}, idx::Int, markov::Int)
    # assumes we have just solved the i'th SP
    @assert idx < M
    sp1 = m.stage_problems[idx, markov]
    for sp2 in m.stage_problems[idx+1,:]
        for key in sp1.ext[:state_vars]
            @assert haskey(sp2.ext[:duals], key)
            chgConstrRHS(sp2.ext[:duals][key], getValue(getVar(sp1, key)))
        end
    end
end

"""
Backward pass for the SDDP algorithm

    backward_pass!(m::SDDPModel)

Inputs:
m        - an SDDPModel

"""
function backward_pass!{M,N,S,T}(m::SDDPModel{M,N,S,T})
    markov = m.init_markov_state==0?transition(m, 1, m.init_markov_state):m.init_markov_state
    old_markov = m.init_markov_state
    old_scenario, scenario=0,0
    for stage=1:M
        # pick a scenario
        scenario=rand(1:S)

        sp = m.stage_problems[stage, markov]      # corresponding stage problem
        load_scenario!(sp, scenario)

        sp.ext[:old_scenario] = (old_markov, old_scenario)  # Lets store the scenario we just came from
        solve!(sp, m)                           # solve
        if stage < M
            old_markov = markov
            old_scenario = scenario
            pass_states_forward!(m, stage, old_markov)# pass forward if necessary
            markov = transition(m, stage, old_markov)           # transition to new scenario

        end
    end
    # Solve all the final stage problems
    for mkv=1:N
        for sc in 1:S
            mkv==markov && sc==scenario && continue

            load_scenario!(m.stage_problems[M,mkv], sc)
            solve!(m.stage_problems[M,mkv], m)
        end
    end

    # Stepping back throught the stages
    for idx=reverse(1:(M-1))
        # Add a cut to that scenario
        add_cut!(m, idx, old_markov)

        # Look up the scenario to step back to
        old_markov, old_scenario = m.stage_problems[idx, old_markov].ext[:old_scenario]

        # Solve all the stage problems
        for sp in m.stage_problems[idx,:]
            for s=1:S
                load_scenario!(sp, s)
                solve!(sp, m)
            end
        end
    end

    # Calculate a new upper bound
    set_valid_bound!(m)
end

function backward_pass!(m::SDDPModel, backward_passes::Int)
    for i=1:backward_passes
        backward_pass!(m)
    end
end

"""
Adds a benders cut to the value to go variable

This function adds a benders cut to the
    add_cut!(m::SDDPModel, stage::Int, scenario::Int)

Inputs:
m        - an SDDPModel
stage    - the stage to add the cut
scenario - the scenario to add the cut to

"""
# function add_cut!{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, markov::Int)
#     # TODO we could tidy this up by collecting terms first
#     #   and then making the cuts
#     sp = m.stage_problems[stage, markov]
#     @defExpr(rhs, sum{
#         get_transition(m, stage, markov, mkv) * (
#             mean(m.stage_problems[stage+1, mkv].ext[:LastObjectives]) +
#             sum{
#                 mean(m.stage_problems[stage+1, mkv].ext[:DualValues][state]) * (
#                     getVar(sp, state) - getRHS(m.stage_problems[stage+1, mkv].ext[:duals][state])
#                 )
#             , state in sp.ext[:state_vars]}
#         )
#     , mkv=1:N}
#     )
#     if m.sense==:Max
#         @addConstraint(sp, sp.ext[:theta] <= rhs)
#     else
#         @addConstraint(sp, sp.ext[:theta] >= rhs)
#     end
# end
function risk_weightings{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, markov::Int)
    P = zeros(S*N)
    i=1
    for mkv=1:N
        for s=1:S
            P[i] = get_transition(m, stage, markov, mkv)/S
            i+=1
        end
    end

    @assert abs(sum(P) - 1) < 1e-5

    w = risk_averse_weightings(vcat([sp.ext[:LastObjectives] for sp in m.stage_problems[stage+1, :]]...), P,  m.beta_quantile, m.sense==:Max)

    Prob = zeros(N,S)
    i=1
    for mkv=1:N
        for s=1:S
            Prob[mkv, s] = m.risk_lambda * P[i] + (1-m.risk_lambda) * w[i]
            i+=1
        end
    end
    return Prob
end

function add_cut!{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, markov::Int)
    sp = m.stage_problems[stage, markov]

    Prob = risk_weightings(m, stage, markov)
    @assert abs(sum(Prob) - 1) < 1e-5
    # if abs(sum(Prob) - 1) > 1e-5
    #     @show stage, markov
    #     @show Prob
    #     @show [sp.ext[:LastObjectives] for sp in m.stage_problems[stage+1, :]]
    #     P = zeros(S*N)
    #     i=1
    #     for mkv=1:N
    #         for s=1:S
    #             P[i] = get_transition(m, stage, markov, mkv)/S
    #             i+=1
    #         end
    #     end
    #
    #     @assert abs(sum(P) - 1) < 1e-5
    #
    #     w = risk_averse_weightings(vcat([sp.ext[:LastObjectives] for sp in m.stage_problems[stage+1, :]]...), P,  m.beta_quantile, m.sense==:Max)
    #     @show w
    #     error()
    # end
    @defExpr(rhs, sum{
        sum{
            Prob[mkv, s] * (
                m.stage_problems[stage+1, mkv].ext[:LastObjectives][s] +
                sum{
                    m.stage_problems[stage+1, mkv].ext[:DualValues][state][s] * (
                        getVar(sp, state) - getRHS(m.stage_problems[stage+1, mkv].ext[:duals][state])
                    )
                , state in sp.ext[:state_vars]}
            )
        ,s=1:S; Prob[mkv, s] > 1e-6}
    , mkv=1:N}
    )

    if m.sense==:Max
        @addConstraint(sp, sp.ext[:theta] <= rhs)
    else
        @addConstraint(sp, sp.ext[:theta] >= rhs)
    end
end

function risk_averse_weightings{T}(x::Vector{T}, p::Vector{T},  ß::Float64=0.5, ismax::Bool=true)
    I = sortperm(x, rev=!ismax)
    y = zeros(length(x))
    q = 0.
    for i in I
        q >=  ß && break
        y[i] = min(p[i],  ß - q)
        q += y[i]
    end
    return y ./ ß
end
risk_averse_weightings{T}(x::Vector{T}, beta::Float64=0.5) = risk_averse_weightings(x, ones(length(x)) / length(x), beta)

function getRHS(c::ConstraintRef{LinearConstraint})
    constr = c.m.linconstr[c.idx]
    sen = JuMP.sense(constr)
    if sen==:range
        error("Range constraints not supported for Scenarios")
    elseif sen == :(==)
        @assert constr.lb == constr.ub
        return constr.ub
    elseif sen == :>=
        return constr.lb
    elseif sen == :<=
        return constr.ub
    end
    error("Sense $(sen) not supported")
end

"""
Solve the model using the SDDP algorithm.
"""
function JuMP.solve{M,N,S,T}(m::SDDPModel{M,N,S,T}; forward_passes=1, backward_passes=1, max_iters=1000, beta_quantile=1, risk_lambda=1)
    print_stats_header()
    # forward_pass!(m, 1)
    # print_stats(m)
    # while not converged
    m.beta_quantile = beta_quantile
    m.risk_lambda = risk_lambda
    i=0
    while i < max_iters
        # Cutting passes
        backward_pass!(m, backward_passes)

        # Simulate
        forward_pass!(m, forward_passes)

        print_stats(m)

        i += backward_passes
    end

end

function print_stats(m::SDDPModel)
    printfmt("({1:8.2f}, {2:8.2f}) | {3:8.2f} | {4:3.2f}\n", m.confidence_interval[1], m.confidence_interval[2], m.valid_bound, 100*rtol(m))
end
function print_stats_header()
    printfmt("{1:22s} | {2:10s} | {3:6s}\n", "Expected Objective", "Valid Bound", "% Gap")
end

end
