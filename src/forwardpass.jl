"""
    forwardpass!(SDDPModel, n)

Perform n forward passes on the SDDPModel
"""
function forwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int, storesamplepoints=false)
    markov = 0                                       # initialise
    setn!(m, n)                                      # store number of forward pass
    for pass = 1:n                                   # for n passes
        markov = m.initial_markov_state              # initial markov state
        if m.initial_markov_state==0
            markov = transition(m, 1, markov)        # transition immediately
        end
        m.forwardstorage.obj[pass] = 0               # reset objective
        for t=1:T
            savemarkov!(m, pass, t, markov)          # store markov state
            sp = subproblem(m, t, markov)            # get subproblem
            load_scenario!(m, sp)                    # realise scenario
            forwardsolve!(sp)                        # solve
            saveobj!(m, pass, getstagevalue(sp))     # store objective
            savex!(m, sp, pass, t)                   # store state
            if storesamplepoints
                addsamplepoint!(m, pass, t, markov)        # store sample points for cutselection
            end
            if t < T
                pass_states!(m, sp, t)               # pass state values forward
                markov = transition(m, t, markov)    # transition
            end
        end
    end
end

# solve the subproblem in a forward pass
function forwardsolve!(sp::Model)
    @assert issubproblem(sp)
    status = solve(sp)
    # Catch case where we aren't optimal
    if status != :Optimal
        sp.internalModelLoaded = false
        status = solve(sp)
        if status != :Optimal
            error("SDDP Subproblems must be feasible. Current status: $(status). I tried rebuilding from the JuMP model but it didn't work...")
        end
    end
end

# load a random scenario
function load_scenario!(m::SDDPModel, sp::Model)
    scenario = sample(m.scenario_probability)
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

function pass_states!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, sp::Model, stage::Int)
    sv = 0.
    for v in stagedata(sp).state_vars
        sv += getvalue(v)
    end
    JuMP.setRHS(stagedata(sp).regularisecons[1], sv)
    JuMP.setRHS(stagedata(sp).regularisecons[2], -sv)
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
