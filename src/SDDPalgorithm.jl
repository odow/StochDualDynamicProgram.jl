"""
This function passes the value of the state variables forward to the next stage

Inputs:
m            - SDDP model object
stage        - the current stage
markov_state - the current markov state
"""
function pass_states_forward!(m::SDDPModel, stage::Int, markov_state::Int)
    # assumes we have just solved the i'th SP

    # Check this is not the last stage
    @assert stage < m.stages

    sp = m.stage_problems[stage, markov_state]

    sv = 0.
    for v in stagedata(sp).state_vars
        sv += getvalue(v)
    end
    JuMP.setRHS(stagedata(sp).regularisecons[1], sv)
    JuMP.setRHS(stagedata(sp).regularisecons[2], -sv)

    # For each of the problems in the next stage
    for next_sp in m.stage_problems[stage+1,:]
        # Check sanity
        @assert length(stagedata(sp).dual_constraints) == length(stagedata(next_sp).dual_constraints) == length(stagedata(next_sp).state_vars)

        # For each state variable
        for i in 1:length(stagedata(next_sp).state_vars)
            JuMP.setRHS(stagedata(next_sp).dual_constraints[i], getvalue(stagedata(sp).state_vars[i]))
        end
    end

    if size(m.stagecuts)[1] > 0
        addsamplepoint!(m.sense, m.stagecuts[stage, markov_state], map(getvalue, stagedata(sp).state_vars))
    end
end

function StatsBase.sample(wv::WeightVec, r::Float64)
    t = r * sum(wv)
    w = values(wv)
    n = length(w)
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end
"""
Backward pass for the SDDP algorithm

    backward_pass!(m::SDDPModel)

Inputs:
m        - an SDDPModel

"""
function load_scenario!(m::SDDPModel, sp::Model)
    scenario = sample(m.scenario_probability)
    load_scenario!(sp, scenario)
    return scenario
end
function load_scenario!(m::SDDPModel, sp::Model, r::Float64)
    scenario = sample(m.scenario_probability, r)
    load_scenario!(sp, scenario)
    return scenario
end

function backward_pass!(m::SDDPModel, cut_selection::Bool, method=NoSelection(), regularisation=NoRegularisation())
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
    for stage=1:(m.stages-1)
        # Get the stage problem
        sp = m.stage_problems[stage, markov]

        # pick a scenario
        scenario = load_scenario!(m, sp)

        # Lets store the scenario we just came from
        stagedata(sp).old_scenario = (old_markov, old_scenario)


        set_regularised_objective!(regularisation, m.sense, sp)
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
    solve_all_stage_problems!(m, m.stages, regularisation)

    # Stepping back throught the stages
    for stage=reverse(1:(m.stages-1))
        # Add a cut to that scenario
        add_cut!(m, stage, old_markov, cut_selection, method)

        # Look up the scenario to step back to
        old_markov, old_scenario = stagedata(m.stage_problems[stage, old_markov]).old_scenario

        # Solve all the stage problems
        solve_all_stage_problems!(m, stage, regularisation)

    end

    # Calculate a new upper bound
    set_valid_bound!(m)
end

"""
This function solves all the stage problems (scenarios x markov states) in a given stage

Inputs:
m     - the SDDP model object
stage - the stage to solve all problems in
"""
function solve_all_stage_problems!(m::SDDPModel, stage::Int, regularisation::Regularisation=NoRegularisation())
    for markov_state in 1:m.markov_states
        set_nonregularised_objective!(regularisation, m.sense, m.stage_problems[stage,markov_state])
        for s=1:m.scenarios
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
function add_cut!(m::SDDPModel, stage::Int, markov_state::Int, cut_selection::Bool=false, method=NoSelection())
    sp = m.stage_problems[stage, markov_state]

    risk_weightings!(m, stage, markov_state)

    @expression(sp, rhs, sum{
        sum{
            m.weightings_matrix[new_markov, scenario] * (
                stagedata(new_sp).objective_value[scenario] +
                sum{
                    stagedata(new_sp).dual_values[i][scenario] * (
                    # (stagedata(new_sp).dual_values[i][scenario] + 1e-10*(rand()-0.5)) * (
                        stagedata(sp).state_vars[i] - getRHS(stagedata(new_sp).dual_constraints[i])
                    )
                , i in 1:length(stagedata(sp).state_vars)}
            )
        ,scenario in 1:m.scenarios; m.weightings_matrix[new_markov, scenario] > 1e-6}
    , (new_markov, new_sp) in enumerate(m.stage_problems[stage+1, :])}
    )

    # for yi in aggregateterms(sp, rhs)
    #     if abs(yi) < 1e-10
    #         warn("Dual Degeneracy detected. When simulating the policy you may find sub-optimal solutions")
    #     end
    # end

    if cut_selection
        if add_cut!(m.sense, sp, rhs, m.stagecuts[stage, markov_state])
            add_cut!(m.sense, sp, rhs, method)
        end
    else
        add_cut!(m.sense, sp, rhs, method)
    end

    m.cuts_filename != nothing && write_cut(m.cuts_filename, sp, stage, markov_state, rhs)

    return
end

function add_cut!(::Type{Val{:Min}}, sp::Model, rhs::JuMP.GenericAffExpr, method::CutSelectionMethod=NoSelection())
    @constraint(sp, stagedata(sp).theta >= rhs)
end
function add_cut!(::Type{Val{:Max}}, sp::Model, rhs::JuMP.GenericAffExpr, method::CutSelectionMethod=NoSelection())
    @constraint(sp, stagedata(sp).theta <= rhs)
end
# function add_cut!(::Type{Val{:Min}}, sp::Model, rhs::JuMP.GenericAffExpr, method::LazyConstraint)
#     @addLazyConstraint(sp, stagedata(sp).theta >= rhs)
# end
# function add_cut!(::Type{Val{:Max}}, sp::Model, rhs::JuMP.GenericAffExpr, method::LazyConstraint)
#     @addLazyConstraint(sp, stagedata(sp).theta <= rhs)
# end
# function corners(cb)
#     x_val = getvalue(x)
#     y_val = getvalue(y)
#     # Allow for some impreciseness in the solution
#     TOL = 1e-6
#     # Check top left, allowing some tolerance
#     if y_val - x_val > 1 + TOL
#         # Cut off this solution
#         println("Solution was in top left, cut it off")
#         # Use the original variables, but not m - cb instead
#         @addLazyConstraint(cb, y - x <= 1)
#     # Check top right
#     elseif y_val + x_val > 3 + TOL
#         # Cut off this solution
#         println("Solution was in top right, cut it off")
#         # Use the original variables, but not m - cb instead
#         @addLazyConstraint(cb, y + x <= 3)
#     end
# end  # End of callback function
#
# # Tell JuMP/Gurobi to use our callback function
# addLazyCallback(m, corners)

function aggregateterms(sp::Model, ex::JuMP.GenericAffExpr)
    data = stagedata(sp)::StageData
    y = zeros(length(data.state_vars))
    for i=1:length(ex.vars)
        for j in 1:length(y)
            if data.state_vars[j].col == ex.vars[i].col
                y[j] += ex.coeffs[i]
            end
        end
    end
    y
end

"""
This function writes the coefficits of a cut to file
"""
function write_cut(filename::ASCIIString, sp::Model, stage::Int, markov_state::Int, rhs::JuMP.GenericAffExpr)
    n = length(rhs.vars)
    open(filename, "a") do f
        write(f, string(stage, ", ", markov_state, ", ", rhs.constant))
        coeffs = aggregateterms(sp, rhs)
        for y in coeffs
            write(f, string(", ", y))
        end
        write(f, "\n")
    end
end

"""
This function recomputes the scenario weightings based on risk aversion parameters

Inputs:
m      - the SDDDP model object
stage  - the current stage
markov - the current markov state
"""
function risk_weightings!(m::SDDPModel, stage::Int, markov::Int)
    # Initialise weightings matrix
    # weightings = zeros(N,S)

    if abs(m.risk_lambda - 1) > 1e-5
        # We are risk averse

        # Initialise storage
        P = zeros(m.scenarios*m.markov_states)

        # Reshape transition matrix appropriately
        i=1
        for mkv=1:m.markov_states
            for s=1:m.scenarios
                P[i] = get_transition(m, stage, markov, mkv)*m.scenario_probability[s]
                i+=1
            end
        end

        # Check sanity
        @assert abs(sum(P) - 1) < 1e-5

        # Construct risk averse weightings for CVar
        w = risk_averse_weightings(vcat([stagedata(sp).objective_value for sp in m.stage_problems[stage+1, :]]...), P,  m.beta_quantile, isa(m.sense, Type{Val{:Max}}))

        # Construct weightings as a convex combination of Expectation (P) and risk averse (w)
        i=1
        for mkv=1:m.markov_states
            for s=1:m.scenarios
                m.weightings_matrix[mkv, s] = m.risk_lambda * P[i] + (1-m.risk_lambda) * w[i]
                i+=1
            end
        end
    else
        # We are just under expectation
        for mkv=1:m.markov_states
            for s=1:m.scenarios
                m.weightings_matrix[mkv, s] = get_transition(m, stage, markov, mkv)*m.scenario_probability[s]
            end
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
function getRHS(c::ConstraintRef)
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

    colval = copy(sp.colVal)

    status = solve(sp)

    # Catch case where we aren't optimal
    if status != :Optimal
        JuMP.writeMPS(sp, "C:/temp/a_infeasible_subproblem_$(myid()).mps")
        sp.internalModelLoaded = false
        status = solve(sp)
        if status != :Optimal
            JuMP.writeMPS(sp, "C:/temp/b_infeasible_subproblem_$(myid()).mps")
            open("C:/temp/b_rhs.csv", "w") do f
                for c in sp.linconstr
                    write(f, "$(string(c)), $(c.lb), $(c.ub)\n)")
                end
            end
            error("SDDP Subproblems must be feasible. Current status: $(status).")
        end
    end

    # Get current scenario
    s = stagedata(sp).current_scenario

    # store the objective value
    stagedata(sp).objective_value[s] = getobjectivevalue(sp)

    # store the dual value for each of the state variables
    for i in 1:length(stagedata(sp).state_vars)
        stagedata(sp).dual_values[i][s] = getdual(stagedata(sp).dual_constraints[i])
    end

    return
end

"""
This function loads a scenario into a Model by changing the RHS

Input:
m        - Stage problem
scenario - the scenario to load
"""
function load_scenario!(sp::Model, scenario::Int)
    # Store old scenario
    stagedata(sp).last_scenario = stagedata(sp).current_scenario

    # for each scenario constraint
    for (constr, Ω) in stagedata(sp).scenario_constraints
        # Look up last scenario.
        if stagedata(sp).last_scenario == 0
            # No old scenario loaded so zero constant
            old_scenario = 0.
        else
            # Look up old constant from scenario set Ω
            old_scenario = Ω[stagedata(sp).last_scenario]
        end

        # Update RHS
        # JuMP.setRHS(constr, getRHS(constr) - old_scenario + Ω[scenario])
        JuMP.setRHS(constr, Ω[scenario])
    end

    # Store new scenario
    stagedata(sp).current_scenario = scenario

    return
end

# Some helper functions
JuMP.getlowerbound(c::ConstraintRef) = c.m.linconstr[c.idx].lb
JuMP.getupperbound(c::ConstraintRef) = c.m.linconstr[c.idx].ub

"""
This function gets the appropriate simulated bound closest to the outer bound
"""
getCloseCIBound(::Type{Val{:Min}}, m::SDDPModel) = m.confidence_interval[1]
getCloseCIBound(::Type{Val{:Max}}, m::SDDPModel) = m.confidence_interval[2]
getCloseCIBound(m::SDDPModel) = getCloseCIBound(m.sense, m)

"""
And vice versa
"""
getFarCIBound(::Type{Val{:Min}}, m::SDDPModel) = m.confidence_interval[2]
getFarCIBound(::Type{Val{:Max}}, m::SDDPModel) = m.confidence_interval[1]
getFarCIBound(m::SDDPModel) = getFarCIBound(m.sense, m)

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
atol(::Type{Val{:Min}}, m::SDDPModel) = getCloseCIBound(m) - getBound(m)
atol(::Type{Val{:Max}}, m::SDDPModel) = getBound(m) - getCloseCIBound(m)
atol(m::SDDPModel) = atol(m.sense, m)

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
function get_true_value(sp::Model)
    @assert is_sp(sp)
    if stagedata(sp).theta != nothing
        return (getobjectivevalue(sp) - getvalue(stagedata(sp).theta))::Float64
    else
        return getobjectivevalue(sp)::Float64
    end
end

"""
Calculate the upper bound of the first stage problem
"""
function set_valid_bound!(m::SDDPModel)
    # Initialise
    obj = 0.0

    if m.init_markov_state == 0
        # Lets average over all first stage probles (markov states x scenarios)
        for i=1:m.markov_states
            for s=1:m.scenarios
                obj += get_transition(m, 1, m.init_markov_state, i) * m.scenario_probability[s] * stagedata(m.stage_problems[1, i]).objective_value[s]::Float64
            end
        end
    else
        # Lets just  average over the scenarios in the initial markov state
        for s=1:m.scenarios
            obj += m.scenario_probability[s] * stagedata(m.stage_problems[1, m.init_markov_state]).objective_value[s]::Float64
        end
    end

    # Update the bound
    decide_set_bound!(m, obj)

    return
end

function decide_set_bound!(::Type{Val{:Min}}, m::SDDPModel, obj)
    if obj > getBound(m)
        setBound!(m, obj)
    end
end
function decide_set_bound!(::Type{Val{:Max}}, m::SDDPModel, obj)
    if obj < getBound(m)
        setBound!(m, obj)
    end
end
decide_set_bound!(m::SDDPModel, obj) = decide_set_bound!(m.sense, m, obj)

function test_and_set_ci!(m::SDDPModel, obj)
    if abs(m.risk_lambda - 1) < 1e-5
        # Not risk averse so Normal Dist CI
        if length(obj) > 1
            setCI!(m, t_test(obj, conf_level=m.QUANTILE))
        else
            setCI!(m, (obj[1], obj[1]))
        end
    else
        # estimate of CVar
        setCI!(m, cvar(obj, m.beta_quantile, m.risk_lambda))
    end
end

@inline store_results!(results::Void, vars, sp, stage, pass, markov, scenario) = nothing
function store_results!(results::Dict{Symbol, Any}, vars, sp, stage, pass, markov, scenario)
    results[:Objective][pass] += get_true_value(sp)         # Add objective
    results[:Current][stage][pass] = get_true_value(sp)
    results[:Scenario][stage][pass] = scenario
    results[:Markov][stage][pass] = markov
    results[:Future][stage][pass] = getobjectivevalue(sp) - get_true_value(sp)
    for v in vars
        results[v][stage][pass] = getvalue(getvariable(sp, v))
    end
end

function forward_pass_kernel!(m::SDDPModel, n::Int, variancereduction::Bool=true, results::Union{Void, Dict{Symbol, Any}}=nothing, vars::Vector{Symbol}=Symbol[])
    # Initialise the objective storage
    obj = zeros(n)

    # For a number of realisations
    pass = 1
    while pass <= n
        Rscenario, Rmarkov = rand(m.stages), rand(m.stages)
        forward_pass_kernel_inner!(m, obj, pass, results, vars, Rscenario, Rmarkov)
        pass += 1
        if variancereduction && ((pass+1) <= n)
            forward_pass_kernel_inner!(m, obj, pass, results, vars, 1 - Rscenario, 1 - Rmarkov)
            pass += 1
        end
    end
    return obj
end

function forward_pass_kernel_inner!(m::SDDPModel, obj::Vector{Float64}, pass::Int, results::Union{Void, Dict{Symbol, Any}}, vars::Vector{Symbol}, Rscenario::Vector{Float64}, Rmarkov::Vector{Float64})
    # Initialise scenario
    markov, old_markov, scenario = 0, 0, 0

    # Transition if need be
    if m.init_markov_state==0
        markov = transition(m, 1, m.init_markov_state, Rmarkov[1])
    else
        markov = m.init_markov_state
    end

    # For all stages
    for stage=1:m.stages
        old_markov = markov

        # realise on scenario
        sp = m.stage_problems[stage, markov]
        scenario = load_scenario!(m, sp, Rscenario[stage])

        # solve
        # set_objective!(Val{:Regularise}, m.sense, sp)
        solve!(sp)

        # Add objective (stage profit only)
        obj[pass] += get_true_value(sp)

        # Save results if necesary
        store_results!(results, vars, sp, stage, pass, markov, scenario)

        # pass forward if necessary
        if stage < m.stages
            pass_states_forward!(m, stage, old_markov)

            # transition to new scenario
            markov = transition(m, stage, old_markov, Rmarkov[stage+1])
        end

    end
end

convergence_pass!(m::SDDPModel, convergence::Convergence) = convergence_pass!(m, convergence.n, convergence.variancereduction)

function convergence_pass!(m::SDDPModel, npasses::Int=1, variancereduction::Bool=true)
    obj = forward_pass_kernel!(m, npasses, variancereduction)

    # set new lower bound
    test_and_set_ci!(m, obj)

    return (rtol(m) < 0., npasses)
end

function convergence_pass!(m::SDDPModel, npasses::Range, variancereduction::Bool)
    OBJ = Float64[]
    for n in npasses
        push!(OBJ, forward_pass_kernel!(m, round(Int, n) - length(OBJ), variancereduction)...)

        test_and_set_ci!(m, OBJ)

        if rtol(m) > 0.
            return (false, n)
        end
        # info("Simulation confidence (n=$n)")
    end
    return (rtol(m) < 0., npasses[end])
end


"""
Simulate SDDP model and return variable solutions contained in [vars]
"""
function serial_simulate(m::SDDPModel, n::Int, vars::Vector{Symbol}=Symbol[]; variancereduction=true)
    results = Dict{Symbol, Any}(:Objective=>zeros(Float64,n))
    for (s, t) in vcat(collect(zip(vars, fill(Any, length(vars)))), [(:Scenario, Int), (:Markov, Int), (:Future, Float64), (:Current, Float64)])
        results[s] = Array(Vector{t}, m.stages)
        for i=1:m.stages
            results[s][i] = Array(t, n)
        end
    end
    forward_pass_kernel!(m, n, variancereduction, results, vars)
    return results
end

"""
Calculate Expected objective
"""
function cvar(x, beta::Float64=1., lambda::Float64=1.)
    @assert beta >= 0 && beta <= 1.
    @assert lambda >= 0 && lambda <= 1.
    cv = lambda * mean(x) + (1 - lambda) * beta_quantile(x, beta)
    return (cv, cv)
end
function beta_quantile(x, beta)
    mean(x[x.<quantile(x, beta)])
end
function t_test(x; conf_level=0.95)
    tstar = quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    SE = std(x)/sqrt(length(x))
    lo, hi = mean(x) + [-1, 1] * tstar * SE
    return (lo, hi)#, mean(x))
end

# @variable(m, x[1:10])
# @variable(m, y)
# @variable(m, z[1:2, [:a, :b]])
#
# results = simulate(m, 1000, x, y, z, z[1,:a] bystage=true)
# results[x][...key...][stage] >> [ ... realisations ...]
# results[y][stage] >> [ ... realisations ...]
# results[z][...key...][stage] >> [ ... realisations ...]
# results[z[1,:a]][stage] >> [ ... realisations ...]
#
# results = simulate(m, 1000, x, y, z, z[1,:a] bystage=false)
# results[x][...key...][...realisation...] >> [ ... stages ...]
# results[y][...realisation...] >> [ ... stages ...]
# results[z][...key...][...realisation...] >> [ ... stages ...]
# results[z[1,:a]][...realisation...] >> [ ... stages ...]
#
# function simulate{M,N,S,T,SEN}(m::SDDPModel{M,N,S,T,SEN}, n::Int, nargs...; bystage=true)
#     results = initialise_results!(Int(M), n, nargs, bystage)
#
#     markov, old_markov = 0, 0
#
#     for pass=1:n
#         if m.init_markov_state==0
#             markov = transition(m, 1, m.init_markov_state)
#         else
#             markov = m.init_markov_state
#         end
#         for stage=1:M
#             old_markov = markov
#
#             scenario = rand(1:S)
#
#             sp = m.stage_problems[stage, markov]  # corresponding stage problem
#             load_scenario!(sp, scenario)
#
#             solve!(sp)                           # solve
#
#             results[:Objective][pass] += get_true_value(sp)         # Add objective
#
#             # pass forward if necessary
#             if stage < M
#                 pass_states_forward!(m, stage, old_markov)
#                 # transition to new scenario
#                 markov = transition(m, stage, old_markov)
#             end
#
#             store_results!(results, sp, stage, pass, nargs, bystage)
#         end
#     end
#     return results
# end
# function initialise_results!(M::Int, n::Int, nargs, bystage)
#     results = Dict{Any, Any}(:Objective=>zeros(Float64,n))
#
#     M = Vector{Float64}[]
#     if bystage
#         M = Vector{Float64}[zeros(M) for i=1:n]
#     else
#         M = Vector{Float64}[zeros(n) for i=1:M]
#     end
#
#     for jump_container in nargs
#         if isa(jump_container, JuMP.Variable)
#             results[jump_container] = copy(M)
#         elseif isa(jump_container, Array{JuMP.Variable})
#             for i in eachindex(jump_container)
#                 results[jump_container][ind2sub(size(jump_container), i)...]=copy(M)
#             end
#         elseif isa(jump_container, JuMP.JuMPArray)
#             for key in keys(jump_container)
#                 results[jump_container][key...] = copy(M)
#             end
#         end
#     end
#
#     return results
# end
#
# function store_results!(results, sp, stage, pass, nargs, bystage)
#     for jump_container in nargs
#         if isa(jump_container, JuMP.Variable)
#             if bystage
#                 results[jump_container][stage][pass] = getvalue(sp, jump_container)
#             else
#                 results[jump_container][pass][stage] = getvalue(sp, jump_container)
#             end
#         elseif isa(jump_container, Array{JuMP.Variable})
#             for i in eachindex(jump_container)
#                 tup = ind2sub(size(jump_container), i)...
#                 if bystage
#                     results[jump_container][tup...][stage][pass] = getvalue(sp, jump_container[tup...])
#                 else
#                     results[jump_container][tup...][pass][stage] = getvalue(sp, jump_container[tup...])
#                 end
#             end
#         elseif isa(jump_container, JuMP.JuMPArray)
#             for key in keys(jump_container)
#                 if bystage
#                     results[jump_container][key...][stage][pass] = getvalue(sp, jump_container[key...])
#                 else
#                     results[jump_container][key...][pass][stage] = getvalue(sp, jump_container[key...])
#                 end
#             end
#         end
#     end
# end
