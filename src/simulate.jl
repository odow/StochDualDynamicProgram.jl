function montecarloestimation{T, M, S, X, TM}(::Type{Val{true}}, m::SDDPModel{T, M, S, X, TM}, n::Int)
    objectives = zeros(n)                                # intial storage for objective
    markov = 0                                           # initialise
    rscenario, rmarkov = zeros(T), zeros(T)              # initialise antithetic storage
    for pass = 1:2:n                                     # for n/2 passes
        @inbounds for i=1:T                              # initalise random variables
            rscenario[i] = rand()                        #    for the scenarios
            rmarkov[i]   = rand()                        #    and markov states
        end
        for antitheticpass = 0:1
            if pass + antitheticpass > n
                break
            end
            markov = m.initial_markov_state              # initial markov state
            if m.initial_markov_state==0                 # no initial state specified
                markov = transition(m, 1, markov,        # transition immediately
                            rmarkov[1])                  #  with the first stored random number
            end
            for t=1:T                                    # for each stage
                sp = subproblem(m, t, markov)            # get subproblem
                load_scenario!(m, sp, rscenario[t])      # realise scenario
                forwardsolve!(sp)                        # solve
                objectives[pass+antitheticpass] += getstagevalue(sp)    # update objective
                if t < T                                 # don't do this for the last stage
                    pass_states!(m, sp, t)               # pass state values forward
                    markov = transition(m, t, markov,    # transition
                                    rmarkov[1+t])        #   since we already used the first
                end
                @inbounds rscenario[t] = 1 - rscenario[t]# get antithetic variate
                @inbounds rmarkov[t]   = 1 - rmarkov[t]  # get antithetic variate
            end
        end
    end
    return objectives
end

function montecarloestimation{T, M, S, X, TM}(::Type{Val{false}}, m::SDDPModel{T, M, S, X, TM}, n::Int)
    objectives = zeros(n)                            # intial storage for objective
    markov = 0                                       # initialise
    for pass = 1:n                                   # for n passes
        markov = m.initial_markov_state              # initial markov state
        if m.initial_markov_state==0                 # no initial state specified
            markov = transition(m, 1, markov)        # transition immediately
        end
        for t=1:T                                    # for each stage
            sp = subproblem(m, t, markov)            # get subproblem
            load_scenario!(m, sp)                    # realise scenario
            forwardsolve!(sp)                        # solve
            objectives[pass] += getstagevalue(sp)    # update objective
            if t < T                                 # don't do this for the last stage
                pass_states!(m, sp, t)               # pass state values forward
                markov = transition(m, t, markov)    # transition
            end
        end
    end
    return objectives
end

function setbound!(log::SolutionLog, m::SDDPModel, bound_convergence)
    old_bound = m.valid_bound
    setbound!(log, m)
    if bound_convergence.after > 0
        if abs(old_bound - m.valid_bound) < bound_convergence.tol
            bound_convergence.n += 1
            if bound_convergence.n > bound_convergence.after
                return false # converged
            end
        else
            bound_convergence.n = 0
        end
    end
    return true # notcongerged
end

"""
    estimatebound!(log, SDDPmodel, convergence, iteration)

Estimate the value of the policy by monte-carlo simulation and using sequential sampling.
"""
function estimatebound!(log::SolutionLog, m::SDDPModel, convergence, iteration::Int, isparallel::Bool, print_level)
    if isparallel
        estimatebound!(log, m, convergence, iteration, parallelmontecarloestimation, print_level)
    else
        estimatebound!(log, m, convergence, iteration, montecarloestimation, print_level)
    end
end

function estimatebound!(log::SolutionLog, m::SDDPModel, convergence, iteration, montecarlofunction::Function, print_level)
    notconverged = true
    ismontecarlo = false
    if convergence.frequency > 0 && mod(iteration, convergence.frequency) == 0
            tic()
            print_level >= PRINTINFO && info("Running out-of-sample Monte Carlo simulation")
            ismontecarlo = true
            obj = copy(getobj(m))
            if length(obj) < convergence.minsamples
                log.simulations += convergence.minsamples - length(obj)
                push!(obj,
                    montecarlofunction(convergence.antithetic, m,
                        convergence.minsamples - length(obj)
                        )...
                )
            end
            ci = estimatebound(obj, convergence.confidencelevel)
            while isconverged(m, ci, m.valid_bound)
                if length(obj) >= convergence.maxsamples
                    if convergence.terminate
                        notconverged = false
                    end
                    break
                end
                oldlength = length(obj)
                push!(obj,
                    montecarlofunction(convergence.antithetic, m,
                        min(
                            convergence.maxsamples-length(obj),
                            convergence.step
                            )
                        )...
                    )
                log.simulations += length(obj) - oldlength
                ci = estimatebound(obj, convergence.confidencelevel)
            end
            setconfidenceinterval!(m, ci)
            log.ci_lower, log.ci_upper = m.confidence_interval
            log.time_forwards += toq()
    else
        setconfidenceinterval!(log, m, 0.95)
    end
    return notconverged, ismontecarlo
end

function initialiseresultsdict(n::Int, T::Int, vars::Vector{Symbol})
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
    return results
end

function serialsimulate{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int, vars::Vector{Symbol})
    results = initialiseresultsdict(n, T, vars)
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

function simulate{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int, vars::Vector{Symbol}=Symbol[]; parallel=false)
    if parallel
        if length(workers()) < 2
            warn("Paralleisation requested but Julia is only running with a single worker. Start julia with `julia -p N` or use the `addprocs(N)` function to load N workers. Running in serial mode.")
            serialsimulate(m, n, vars)
        else
            parallelsimulate(m, n, vars)
        end
    else
        serialsimulate(m, n, vars)
    end
end

function historicalsimulation{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, vars::Vector{Symbol}=Symbol[]; markov=ones(Int, m.stages), kwargs...)
    n, scenario, obj = 1, 0, 0.          # Initialise storage
    results = initialiseresultsdict(n, T, vars)
    for t=1:T                            # For all stages
        sp = subproblem(m, t, markov[t]) # get subproblem
        s = load_scenario!(m, sp)        # realise scenario
        data = stagedata(sp)
        for (key, series) in kwargs
            @assert haskey(data.scenario_constraint_names, key)
            cidx = data.scenario_constraint_names[key]
            JuMP.setRHS(data.scenario_constraints[cidx][1], series[t])
        end
        forwardsolve!(sp)                # solve
        store_results!(results, vars, sp, t, 1, markov[t], s)
        t < T && pass_states!(m, sp, t)  # pass state values forward
    end
    results
end
