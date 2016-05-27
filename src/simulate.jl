function simulate{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int, vars::Vector{Symbol}=Symbol[]; variancereduction=true)
    results = Dict{Symbol, Any}(:Objective=>zeros(Float64,n))
    for (s, ty) in vcat(
                    collect(zip(vars, fill(Any, length(vars)))),
                    [(:Scenario, Int), (:Markov, Int), (:Future, Float64), (:Current, Float64)]
                    )
        results[s] = Array(Vector{ty}, T)
        for t=1:T
            results[s][t] = Array(ty, n)
        end
    end
    i = 0                                  # initialise
    for pass = 1:n                         # for n passes
        i = m.initial_markov_state         # initial markov state
        if m.initial_markov_state==0
            i = transition(m, 1, i)        # transition immediately
        end
        for t=1:T
            sp = subproblem(m, t, i)       # get subproblem
            s = load_scenario!(m, sp)      # realise scenario
            forwardsolve!(sp)              # solve
            store_results!(results, vars, sp, t, pass, i, s)
            if t < T
                pass_states!(m, sp, t)     # pass state values forward
                i = transition(m, t, i)    # transition
            end
        end
    end
    return results
end

@inline store_results!(results::Void, vars, sp, stage, pass, markov, scenario) = nothing
function store_results!(results::Dict{Symbol, Any}, vars, sp, stage, pass, markov, scenario)
    results[:Objective][pass] += getstagevalue(sp)
    results[:Current][stage][pass] = getstagevalue(sp)
    results[:Scenario][stage][pass] = scenario
    results[:Markov][stage][pass] = markov
    results[:Future][stage][pass] = getobjectivevalue(sp) - getstagevalue(sp)
    for v in vars
        results[v][stage][pass] = getvalue(getvariable(sp, v))
    end
end
