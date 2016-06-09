"""
    forwardpass!(SDDPModel, n)

Perform n forward passes on the SDDPModel
"""
function forwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int, cutselection::CutSelectionMethod, forwardpass::ForwardPass)
    markov = 0                                       # initialise
    for pass = 1:n                                   # for n passes
        markov = m.initial_markov_state              # initial markov state
        if m.initial_markov_state==0                 # no initial state specified
            markov = transition(m, 1, markov)        # transition immediately
        end
        m.forwardstorage.obj[pass] = 0               # reset objective
        for t=1:T                                    # for each stage
            savemarkov!(m, pass, t, markov)          # store markov state
            sp = subproblem(m, t, markov)            # get subproblem
            load_scenario!(m, sp)                    # realise scenario
            regularisedsolve!(X, sp, forwardpass.regularisation)
            saveobj!(m, pass, getstagevalue(sp))     # store objective
            savex!(m, sp, pass, t)                   # store state
            addsamplepoint!(m, cutselection,         # store sample points for
                    pass, t, markov)                 #    cutselection
            if t < T                                 # don't do this for the last stage
                pass_states!(m, sp, t, forwardpass.regularisation) # pass state values forward
                markov = transition(m, t, markov, forwardpass.uniformsampling)    # transition
            end
        end
    end
end

# a wrapper for forward pass timings
function forwardpass!(log::SolutionLog, m::SDDPModel, iterations, cutselection, forward_pass, isparallel)
    tic()
    # number of scenarios to sample on forward pass
    nscenarios = getscenarios(forward_pass.scenarios, iterations)
    if isparallel
        parallelforwardpass!(m, nscenarios, cutselection, forward_pass)
    else
        resizeforwardstorage!(m, nscenarios) # resize storage for forward pass
        forwardpass!(m, nscenarios, cutselection, forward_pass)
    end
    log.simulations += nscenarios
    log.time_forwards += toq()
    return nscenarios
end

"""
    getscenarios(forwardscenarios, iteration)

This function gets the number of scenarios to sample in iteration `iteration`.
"""
getscenarios(forwardscenarios::Int, iteration::Int) = forwardscenarios
function getscenarios(forwardscenarios::AbstractArray{Int,1}, iteration::Int)
    if iteration < length(forwardscenarios)
        return forwardscenarios[iteration]
    else
        return forwardscenarios[end]
    end
end

function load_scenario!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, sp::Model, importance::Bool)
    if importance
        load_scenario!(m, sp, rand(1:S))
    else
        return load_scenario!(m, sp)
    end
end

function transition{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t, markov, importance::Bool)
    if importance
        return rand(1:M)
    else
        return transition(m, t, markov)
    end
end

function regularisedsolve!(sense::Sense, sp::Model, regularisation::Regularisation)
    set_regularised_objective!(regularisation, sense, sp)
    forwardsolve!(sp)                        # solve
    set_nonregularised_objective!(regularisation, sense, sp)
end
regularisedsolve!(sense::Sense, sp::Model, regularisation::NoRegularisation) = forwardsolve!(sp)

addsamplepoint!(m::SDDPModel, cs::LevelOne, pass, t, markov) = addsamplepoint!(m, pass, t, markov)
addsamplepoint!(m::SDDPModel, cs::CutSelectionMethod, pass, t, markov) = nothing

function resizeforwardstorage!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int)
    resize!(m.forwardstorage.obj, n)
    nx = length(stagedata(m, 1,1).state_vars)
    for i=(m.forwardstorage.n+1):n
        push!(m.forwardstorage.x, zeros(nx, T))
        push!(m.forwardstorage.W, zeros(Int, T))
    end
    setn!(m, n)
end

function force_resizeforwardstorage!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int)
    m.forwardstorage = ForwardPassData(n)
    nx = length(stagedata(m, 1,1).state_vars)
    for i=1:n
        m.forwardstorage.x[i] = zeros(Float64, (nx, T))
        m.forwardstorage.W[i] = zeros(Int, T)
    end
end

# solve the subproblem in a forward pass
function forwardsolve!(sp::Model)
    @assert issubproblem(sp)
    status = solve(sp)
    # Catch case where we aren't optimal
    if status != :Optimal
        warn("SDDP subproblem not optimal (stats=$(status)). Assuming numerical infeasibility so rebuilding model from stored cuts.")
        sp.internalModelLoaded = false
        status = solve(sp)
        if status != :Optimal
            JuMP.writeMPS(sp, "subproblem_proc$(myid())_$(randstring(8)).mps")
            error("SDDP Subproblems must be feasible. Current status: $(status). I tried rebuilding from the JuMP model but it didn't work so I wrote you an MPS file.")
        end
    end
end

# load a random scenario
function load_scenario!(m::SDDPModel, sp::Model)
    scenario = sample(m.scenario_probability)
    load_scenario!(sp, scenario)
    return scenario
end

function load_scenario!(m::SDDPModel, sp::Model, rnd::Float64)
    scenario = sample(m.scenario_probability, rnd)
    load_scenario!(sp, scenario)
    return scenario
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

# Load a specific scenario
function load_scenario!(sp::Model, scenario::Int)
    for (constr, Ω) in stagedata(sp).scenario_constraints
        JuMP.setRHS(constr, Ω[scenario])
    end
end

function setregularisation!(::LinearRegularisation, sp::Model)
    sv = 0.
    for v in stagedata(sp).state_vars
        sv += getvalue(v)
    end
    JuMP.setRHS(stagedata(sp).regularisecons[1], sv)
    JuMP.setRHS(stagedata(sp).regularisecons[2], -sv)
end
setregularisation!(::Regularisation, sp::Model) = nothing

function pass_states!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, sp::Model, stage::Int, regularisation::Regularisation=NoRegularisation())
    setregularisation!(regularisation, sp)

    # For each of the problems in the next stage
    for next_sp in subproblems(m, stage+1)
        # Check sanity
        @assert length(stagedata(sp).dual_constraints) == length(stagedata(next_sp).dual_constraints) == length(stagedata(next_sp).state_vars)
        # For each state variable
        for i in 1:length(stagedata(next_sp).state_vars)
            JuMP.setRHS(stagedata(next_sp).dual_constraints[i], getvalue(stagedata(sp).state_vars[i]))
        end
    end
end
