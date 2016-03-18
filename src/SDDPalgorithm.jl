"""
This function passes the value of the state variables forward to the next stage

Inputs:
m            - SDDP model object
stage        - the current stage
markov_state - the current markov state
"""
function pass_states_forward!{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, markov_state::Int)
    # assumes we have just solved the i'th SP

    # Check this is not the last stage
    @assert stage < M

    sp = m.stage_problems[stage, markov_state]

    # For each of the problems in the next stage
    for next_sp in m.stage_problems[stage+1,:]
        # For each state variable
        for state in next_sp.ext[:state_vars]
            # Check sanity
            @assert haskey(next_sp.ext[:dual_constraints], state)

            # Change RHS on dummy constraint
            chgConstrRHS(next_sp.ext[:dual_constraints][state], getValue(getVar(sp, state)))
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
    # Initialise markov state
    markov = 0
    if m.init_markov_state==0
        markov = transition(m, 1, m.init_markov_state)
    else
        markov = m.init_markov_state
    end

    # Initialise some storage
    old_markov, old_scenario, scenario=m.init_markov_state, 0,0

    # For all stages going forwards
    for stage=1:(M-1)
        # Get the stage problem
        sp = m.stage_problems[stage, markov]

        # pick a scenario
        scenario=rand(1:S)
        load_scenario!(sp, scenario)

        # Lets store the scenario we just came from
        sp.ext[:old_scenario] = (old_markov, old_scenario)

        # solve
        solve!(sp)

        # Keep track of old markov state and scenario
        old_markov, old_scenario = markov, scenario

        # pass state variables forward
        pass_states_forward!(m, stage, old_markov)

        # transition to new markov state
        markov = transition(m, stage, old_markov)
    end

    # Solve all the final stage problems
    solve_all_stage_problems!(m, M)

    # Stepping back throught the stages
    for stage=reverse(1:(M-1))
        # Add a cut to that scenario
        add_cut!(m, stage, old_markov)

        # Look up the scenario to step back to
        old_markov, old_scenario = m.stage_problems[stage, old_markov].ext[:old_scenario]

        # Solve all the stage problems
        solve_all_stage_problems!(m, stage)
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
This function solves all the stage problems (scenarios x markov states) in a given stage

Inputs:
m     - the SDDP model object
stage - the stage to solve all problems in
"""
function solve_all_stage_problems!{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int)
    for markov_state in 1:N
        for s=1:S
            load_scenario!(m.stage_problems[stage,markov_state], s)
            solve!(m.stage_problems[stage,markov_state])
        end
    end
end

"""
Adds a benders cut to the value to go variable

This function adds a benders cut to the
    add_cut!(m::SDDPModel, stage::Int, markov_state::Int)

Inputs:
m            - an SDDPModel object
stage        - the stage to add the cut
markov_state - the scenario to add the cut to
"""
function add_cut!{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, markov_state::Int)
    sp = m.stage_problems[stage, markov_state]

    risk_weightings!(m, stage, markov_state)

    @defExpr(rhs, sum{
        sum{
            m.weightings_matrix[new_markov, scenario] * (
                new_sp.ext[:objective_value][scenario] +
                sum{
                    new_sp.ext[:dual_values][state][scenario] * (
                        getVar(sp, state) - getRHS(new_sp.ext[:dual_constraints][state])
                    )
                , state in sp.ext[:state_vars]}
            )
        ,scenario in 1:S; m.weightings_matrix[new_markov, scenario] > 1e-6}
    , (new_markov, new_sp) in enumerate(m.stage_problems[stage+1, :])}
    )

    if m.sense==:Max
        @addConstraint(sp, sp.ext[:theta] <= rhs)
    else
        @addConstraint(sp, sp.ext[:theta] >= rhs)
    end
end

"""
This function recomputes the scenario weightings based on risk aversion parameters

Inputs:
m      - the SDDDP model object
stage  - the current stage
markov - the current markov state
"""
function risk_weightings!{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, markov::Int)
    # Initialise weightings matrix
    # weightings = zeros(N,S)

    if abs(m.risk_lambda - 1) > 1e-5
        # We are risk averse

        # Initialise storage
        P = zeros(S*N)

        # Reshape transition matrix appropriately
        i=1
        for mkv=1:N
            for s=1:S
                P[i] = get_transition(m, stage, markov, mkv)/S
                i+=1
            end
        end

        # Check sanity
        @assert abs(sum(P) - 1) < 1e-5

        # Construct risk averse weightings for CVar
        w = risk_averse_weightings(vcat([sp.ext[:objective_value] for sp in m.stage_problems[stage+1, :]]...), P,  m.beta_quantile, m.sense==:Max)

        # Construct weightings as a convex combination of Expectation (P) and risk averse (w)
        i=1
        for mkv=1:N
            for s=1:S
                m.weightings_matrix[mkv, s] = m.risk_lambda * P[i] + (1-m.risk_lambda) * w[i]
                i+=1
            end
        end
    else
        # We are just under expectation
        for mkv=1:N
            m.weightings_matrix[mkv, :] = get_transition(m, stage, markov, mkv)/S
        end
    end

    # Check sanity
    @assert abs(sum(m.weightings_matrix) - 1) < 1e-5

    return
end

"""
This function computes the risk averse weightings for CVar given the beta quantile

Inputs:
x     - objective values for scenarios
p     - probability density for scenarios
ß     - CVar beta quantile
ismax - (true/false) problem is maximisation
"""
function risk_averse_weightings{T}(x::Vector{T}, p::Vector{T},  ß::Float64=0.5, ismax::Bool=true)
    # sort scenarios
    I = sortperm(x, rev=!ismax)

    # Initialise output
    y = zeros(length(x))

    # Quantile collected so far
    q = 0.

    # For each scenario in order
    for i in I
        # If we have collected more than beta quantile stop
        q >=  ß && break

        # Else lets take the biggest proportion of the scenario possible
        y[i] = min(p[i],  ß - q)

        # Update total quantile
        q += y[i]
    end

    # Normalise weights
    return y ./ ß
end
# If we don't give a probability support vector assume uniform
risk_averse_weightings{T}(x::Vector{T}, beta::Float64=0.5) = risk_averse_weightings(x, ones(length(x)) / length(x), beta)

"""
Get the RHS of a linear JuMP constraint
"""
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
Solve a StageProblem.

    solve!(sp::Model)
"""
function solve!(sp::Model)
    @assert is_sp(sp)

    status = solve(sp)

    # Catch case where we aren't optimal
    if status != :Optimal
        error("SDDP Subproblems must be feasible. Current status: $(status).")
    end

    # Get current scenario
    s = sp.ext[:CurrentScenario]

    # store the objective value
    sp.ext[:objective_value][s] = getObjectiveValue(sp)

    # store the dual value for each of the state variables
    for v in sp.ext[:state_vars]
        sp.ext[:dual_values][v][s] = getDual(sp.ext[:dual_constraints][v])
    end

    return
end

"""
This function loads a scenario into a Model by changing the RHS

Input:
m        - Stage problem
scenario - the scenario to load
"""
function load_scenario!(m::Model, scenario::Int)
    # Store old scenario
    m.ext[:LastScenario] = m.ext[:CurrentScenario]

    # for each scenario constraint
    for (constr, Ω) in m.ext[:scenario_constraints]
        # Look up last scenario.
        if m.ext[:LastScenario] == 0
            # No old scenario loaded so zero constant
            old_scenario = 0.
        else
            # Look up old constant from scenario set Ω
            old_scenario = Ω[m.ext[:LastScenario]]
        end

        # Update RHS
        chgConstrRHS(constr, getRHS(constr) - old_scenario + Ω[scenario])
    end

    # Store new scenario
    m.ext[:CurrentScenario] = scenario

    return
end

# Some helper functions
JuMP.getLower(c::ConstraintRef) = c.m.linconstr[c.idx].lb
JuMP.getUpper(c::ConstraintRef) = c.m.linconstr[c.idx].ub

"""
This function gets the appropriate simulated bound closest to the outer bound
"""
function getCloseCIBound(m::SDDPModel)
    if m.sense==:Max
        m.confidence_interval[2]
    else
        m.confidence_interval[1]
    end
end
"""
And vice versa
"""
function getFarCIBound(m::SDDPModel)
    if m.sense==:Max
        m.confidence_interval[1]
    else
        m.confidence_interval[2]
    end
end
"""
Get the outer bound for the model
"""
getBound(m::SDDPModel) = m.valid_bound

"""
Set the confidence interval
"""
setCI!{T<:Real}(m::SDDPModel, v::Tuple{T,T}) = (m.confidence_interval = v)
"""
Set the outer bound
"""
setBound!{T<:Real}(m::SDDPModel, v::T) = (m.valid_bound = v)

"""
Actual tolerance of the solution

Defined as [Outer bound - closest simulated bound]
"""
function atol(m)
    if m.sense==:Max
        getBound(m) - getCloseCIBound(m)
    else
        getCloseCIBound(m) - getBound(m)
    end
end

"""
Relative tolerance of the solution

Defined as [Outer bound - closest simulated bound] / [Outer bound]
"""
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
    # Initialise
    obj = 0.0

    if m.init_markov_state == 0
        # Lets average over all first stage probles (markov states x scenarios)
        for i=1:N
            for s=1:S
                obj += get_transition(m, 1, m.init_markov_state, i) * m.stage_problems[1, i].ext[:objective_value][s]
            end
        end
        obj /= S
    else
        # Lets just  average over the scenarios in the initial markov state
        for s=1:S
            obj += m.stage_problems[1, m.init_markov_state].ext[:objective_value][s]
        end
        obj /= S
    end

    # Update the bound
    if (m.sense==:Max && obj < getBound(m)) || (m.sense==:Min && obj > getBound(m))
        setBound!(m, obj)
    end
end

"""
Forward pass of the SDDP algorithm
"""
function forward_pass!{M,N,S,T}(m::SDDPModel{M,N,S,T}, npasses::Int=1)
    # Initialise scenario
    scenario = 0

    # Initialise the objective storage
    obj = zeros(npasses)

    # For a number of realisations
    for pass=1:npasses
        # Transition if need be
        if m.init_markov_state==0
            markov = transition(m, 1, m.init_markov_state)
        else
            markov = m.init_markov_state
        end

        # For all stages
        for stage=1:M
            old_markov = markov

            # realise on scenario
            sp = m.stage_problems[stage, markov]
            scenario = rand(1:S)
            load_scenario!(sp, scenario)

            # solve
            solve!(sp)

            # Add objective (stage profit only)
            obj[pass] += get_true_value(sp)

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
        # Not risk averse so Normal Dist CI
        if npasses > 1
            _obj = t_test(obj, conf_level=m.QUANTILE)
        else
            _obj = (obj[1], obj[1])
        end
    else
        # Bootstrap estimate of CVar
        _obj = cvar_bootstrap(obj, m.beta_quantile, m.risk_lambda)
    end

    # Update lower bound
    if (m.sense==:Max && _obj[1] > getFarCIBound(m)) || (m.sense==:Min && _obj[2] < getFarCIBound(m))
        setCI!(m, _obj)
    end

end

function cvar_bootstrap(x, beta=1., lambda=1., n=1000)
    y = zeros(100)#div(length(x), 100))
    z = zeros(n)
    for i=1:n
        rand!(y, x)
        z[i] = cvar(y, beta, lambda)
    end
    t_test(z)
end

"""
Calculate Expected objective
"""
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

            solve!(sp)                           # solve

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
