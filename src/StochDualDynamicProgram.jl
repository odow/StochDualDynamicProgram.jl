module StochDualDynamicProgram

importall JuMP
using MathProgBase
using Clp
using Formatting
using Distributions

export SDDPModel,
    @defStateVar, @defValueToGo,
    addStageProblem!

type SDDPModel{M,N,T}
    stage_problems::Array{JuMP.Model, 2}
    transition::Array{Float64, T}
    initial_scenario::Int
    lower_bound::Tuple{Float64, Float64}
    upper_bound::Float64
    status::Symbol
    QUANTILE::Float64
    ATOL::Float64
    RELTOL::Float64
    LPSOLVER::MathProgBase.AbstractMathProgSolver
end

"""
Creates an SDDPModel type

    SDDPModel([;stages, scenarios, transition, initial_scenario, conf_level, lpsolver, mipsolver, abs_tol, rel_tol])

Inputs:
stages     - the number of stages
scenarios  - the number of scenarios
transition - the transition matrix. Can be N*N where N is the number of scenarios, or M*N*N where M is the number of stages.
initial_scenario - the scenario at time 0. Model transitions at the start of time period 1
conf_level - confidence level for the lower bound
lpsolver   - AbstractMathProgBase solver capable of returning dual variables. Defaults to Clp.
abs_tol    - Absolute tolerance convergence criteria. DDP halts when ub - lb < abs_tol
rel_tol    - Relative tolerance convergence criteria. DDP halts when (ub - lb) / lb < rel_tol
"""
function SDDPModel(;
    stages=0,
    scenarios=1,
    transition=nothing,
    initial_scenario=0,
    conf_level=0.95,
    lpsolver=ClpSolver(),
    abs_tol=1e-8,
    rel_tol=1e-8)
    if transition==nothing
        # Uniform transition matrix
        transition = 1 / scenarios * ones(scenarios, scenarios)
    end
    T = length(size(transition))
    # Check valid transition matrix

    if T==2
        @assert size(transition)[1] == size(transition)[2] == scenarios
        for i=1:scenarios
            @assert abs(sum(transition[i,:]) - 1.) < 1e-10
        end
    elseif T==3
        @assert size(transition)[1] == stages
        @assert size(transition)[2] == size(transition)[3] == scenarios
        for i=1:size(transition)[1]
            for j=1:scenarios
                @assert abs(sum(transition[i,j,:]) - 1.) < 1e-10
            end
        end
    end
    SDDPModel{stages,scenarios,T}(Array(JuMP.Model, (stages, scenarios)), transition, initial_scenario, (-Inf, -Inf), Inf, :Unconverged, conf_level, abs_tol, rel_tol, lpsolver)
end
function Base.copy{M,N,T}(m::SDDPModel{M,N,T})
    SDDPModel{M,N,T}(deepcopy(m.stage_problems), copy(m.transition), m.initial_scenario, m.lower_bound, m.upper_bound, m.status, m.QUANTILE, m.ATOL, m.RELTOL, m.LPSOLVER)
end

function transition{M,N,T}(m::SDDPModel{M,N,T}, stage::Int, scenario::Int)
    r = rand()
    for i=1:N
        k = get_transition(m, stage, scenario, i)
        if r <= k
            return i
        else
            r -= k
        end
    end
    return N
end
# transition{M,N}(m::SDDPModel{M,N,2}, scenario::Int) = transition(m, 1, scenario)

function get_transition{M,N}(m::SDDPModel{M,N,2}, stage::Int, s1::Int, s2::Int)
    s1==0?1/N : m.transition[s1,s2]
end
# get_transition{M,N}(m::SDDPModel{M,N,2}, s1::Int, s2::Int) = get_transition(m, 1, s1, s2)
function get_transition{M,N}(m::SDDPModel{M,N,3}, stage::Int, s1::Int, s2::Int)
    s1==0?1/N : m.transition[stage, s1,s2]
end

"""
Instaniates a new StageProblem which is a JuMP.Model object with an extension
dictionary.
"""
function StageProblem()
    sp = Model()
    sp.ext[:is_sp] = true
    sp.ext[:state_vars] = Symbol[]
    sp.ext[:duals] = Dict{Symbol, ConstraintRef}()
    sp.ext[:theta] = nothing
    sp.ext[:old_scenario] = 0
    return sp
end

"""
Add a new stage problem to the SDDPModel.

Usage

    addStageProblem!(m, stage, scenario) do sp
        @defVar(sp, x)
    end

"""
function addStageProblem!{M,N,T}(build_model!::Function, m::SDDPModel{M,N,T}, stage_idx::Int, scenario_idx::Int)
    sp = StageProblem()
    build_model!(sp)

    # remove last theta
    if stage_idx == M && sp.ext[:theta] != nothing
        @addConstraint(sp, sp.ext[:theta] == 0)
    end

    m.stage_problems[stage_idx, scenario_idx] = sp
end
addStageProblem!(build_model!::Function, m::SDDPModel, stage_idx::Int) = addStageProblem!(build_model!, m, stage_idx, 1)

# Check this is a StageProblem
is_sp(m::JuMP.Model) = haskey(m.ext, :is_sp) && m.ext[:is_sp]
is_sp(m) = false

"""
A helper function for the JuMP.getValue so it can extract the values from and SDDP problem.
"""
function JuMP.getValue(m::SDDPModel, v::Symbol)
    x = Array(Any, length(m.stage_problems))
    for (idx, sp) in enumerate(m.stage_problems)
        try
            x[idx] = getValue(getVar(sp, v))
        catch
        end
    end
    if issubtype(typeof(x[1]), Number)
        return convert(Vector{Float64}, x)
    else
        # Deal with more complicated values
        return x
    end
end

# Some helper functions
JuMP.getLower(c::ConstraintRef) = c.m.linconstr[c.idx].lb
JuMP.getUpper(c::ConstraintRef) = c.m.linconstr[c.idx].ub
JuMP.getLower(m::SDDPModel) = m.lower_bound[2]
getLowerBnd(m::SDDPModel) = m.lower_bound[1]
JuMP.getUpper(m::SDDPModel) = m.upper_bound
setLower!{T<:Real}(m::SDDPModel, v::Tuple{T,T}) = (m.lower_bound = v)
setUpper!{T<:Real}(m::SDDPModel, v::T) = (m.upper_bound = v)
atol(m) = getUpper(m) - getLower(m)
rtol(m) = abs(getUpper(m)) != Inf?abs(atol(m) / getUpper(m)):Inf

"""
Get the value of the current stage (i.e. not including the value to go)
"""
function get_true_value(m::Model)
    @assert is_sp(m)
    if m.ext[:theta] != nothing
        return getObjectiveValue(m) - getValue(m.ext[:theta])
    else
        return getObjectiveValue(m)
    end
end

"""
Calculate the upper bound of the first stage problem
"""
function set_upper_bound!{M,N,T}(m::SDDPModel{M,N,T})
    obj = 0.0
    if m.initial_scenario == 0
        for i=1:N
            obj += get_transition(m, 1, m.initial_scenario, i) * getObjectiveValue(m.stage_problems[1, i])
        end
    else
        obj += getObjectiveValue(m.stage_problems[1, m.initial_scenario])
    end
    if obj < getUpper(m)
        setUpper!(m, obj)
    end
end

"""
Solve a StageProblem.

    solve!(sp::Model[, relaxation])
"""
function solve!(sp::Model, m::SDDPModel)
    # Catch case where we aren't optimal
    if solve(sp) != :Optimal
        println("Model:")
        display(sp)
        println("\nConstraints:")
        display(sp.linconstr)
        println()
        error("Unable to continue. Not solved to optimality")
    end
end

"""
Forward pass of the MISDDP algorithm
"""
function forward_pass!{M,N,T}(m::SDDPModel{M,N,T}, npasses::Int=1, print_info::Bool=false)
    scenario = 0                                    # Initialise variable
    obj = zeros(npasses)                            # Initialise the objective values

    for pass=1:npasses
        if m.initial_scenario==0
            scenario = transition(m, 1, m.initial_scenario)
        else
            scenario = m.initial_scenario
        end
        for stage=1:M

            old_scenario = scenario

            sp = m.stage_problems[stage, scenario]  # corresponding stage problem
            sp.ext[:old_scenario] = old_scenario    # Lets store the scenario we just came from
            solve!(sp, m)                           # solve
            obj[pass] += get_true_value(sp)         # Add objective

            # pass forward if necessary
            if stage < M
                pass_states_forward!(m, stage, scenario)

                # transition to new scenario
                scenario = transition(m, stage, scenario)
            end
        end
    end

    # set new lower bound
    if npasses > 1
        _obj = t_test(obj, conf_level=m.QUANTILE)
    else
        _obj = (obj[1], obj[1])
    end
    if _obj[1] > getLowerBnd(m)
        setLower!(m, _obj)
    end

end

function t_test(x::Vector; conf_level=0.95)
    tstar = 1.96 # quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    SE = std(x)/sqrt(length(x))
    lo, hi = mean(x) + [-1, 1] * tstar * SE
    return (lo, hi)#, mean(x))
end

function pass_states_forward!{M,N,T}(m::SDDPModel{M,N,T}, idx::Int, scenario::Int)
    # assumes we have just solved the i'th SP
    @assert idx < M
    sp1 = m.stage_problems[idx, scenario]
    for sp2 in m.stage_problems[idx+1,:]
        for key in sp1.ext[:state_vars]
            @assert haskey(sp2.ext[:duals], key)
            chgConstrRHS(sp2.ext[:duals][key], getValue(getVar(sp1, key)))
        end
    end
end

"""
Backward pass for the MISDDP algorithm

    backward_pass!(m::SDDPModel)

Inputs:
m        - an SDDPModel

"""
function backward_pass!{M,N,T}(m::SDDPModel{M,N,T})
    scenario = m.initial_scenario==0?transition(m, 1, m.initial_scenario):m.initial_scenario
    old_scenario = m.initial_scenario
    for stage=1:M
        sp = m.stage_problems[stage, scenario]  # corresponding stage problem
        sp.ext[:old_scenario] = old_scenario    # Lets store the scenario we just came from
        solve!(sp, m)                           # solve
        if stage < M
            old_scenario = scenario
            pass_states_forward!(m, stage, old_scenario) # pass forward if necessary
            scenario = transition(m, stage, old_scenario)      # transition to new scenario
        end
    end
    # Solve all the final stage problems
    for sp=1:N
        sp==scenario && continue
        solve!(m.stage_problems[M,sp], m)
    end

    # Stepping back throught the stages
    for idx=reverse(1:(M-1))
        # Add a cut to that scenario
        add_cut!(m, idx, old_scenario)

        # Look up the scenario to step back to
        old_scenario = m.stage_problems[idx, old_scenario].ext[:old_scenario]

        # Solve all the stage problems
        for sp=1:N
            solve!(m.stage_problems[idx,sp], m)
        end
    end

    # Calculate a new upper bound
    set_upper_bound!(m)
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
function add_cut!(m::SDDPModel, stage::Int, scenario::Int)
    # Get the stage problem to add the cut to
    sp = m.stage_problems[stage, scenario]

    @addConstraint(sp,
        sp.ext[:theta] <= sum{
            get_transition(m, stage, scenario, i) * (
                getObjectiveValue(s2) +
                sum{getDual(s2.ext[:duals][state]) * (
                    getVar(sp, state) - getLower(s2.ext[:duals][state])
                    ), state=sp.ext[:state_vars]
                }
            ), (i, s2)=enumerate(m.stage_problems[stage+1,:])
        }
    )
end

"""
Solve the model using the SDDP algorithm.

Parameters:
- verbosity
    0 = nothing
    1 = print iteration status
    2 = info messages

Returns:
- m::SDDPModel
    the converged SDDPModel
"""
function JuMP.solve{M,N,T}(m::SDDPModel{M,N,T}; verbosity=1, forward_passes=1, backward_passes=1)
    print_stats_header()
    forward_pass!(m, 1)
    print_stats(m)
    # while not converged
    while atol(m) > m.ATOL && rtol(m) > m.RELTOL
        # Cutting passes
        backward_pass!(m, backward_passes)

        # Simulate
        forward_pass!(m, forward_passes)

        print_stats(m)
    end

end

function print_stats(m::SDDPModel)
    printfmt("({1:8.2f}, {2:8.2f}) | {3:8.2f} | {4:3.2f}\n", m.lower_bound[1], m.lower_bound[2], m.upper_bound, 100*rtol(m))
end
function print_stats_header()
    printfmt("{1:22s} | {2:10s} | {3:6s}\n", "LB", "UB", "% Gap")
end


end
