function simulate(m::SDDPModel, n::Int, vars::Vector{Symbol}=[]; parallel=false)
    if parallel && length(workers()) < 2
        warn("Paralleisation requested but Julia is only running with a single processor. Start julia with `julia -p N` or use the `addprocs(N)` function to load N workers. Running in serial mode.")
        parallel = false
    end
    if parallel
        results = parallel_simulate(m, n, vars)
    else
        results = serial_simulate(m, n, vars)
    end
    results
end

function simulate(m::SDDPModel, vars::Vector{Symbol}=Symbol[]; markov=ones(Int, m.stages), kwargs...)
    n=1
    results = Dict{Symbol, Any}(:Objective=>zeros(Float64,n))
    for (s, t) in vcat(collect(zip(vars, fill(Any, length(vars)))), [(:Scenario, Int), (:Markov, Int), (:Future, Float64), (:Current, Float64)])
        results[s] = Array(Vector{t}, m.stages)
        for i=1:m.stages
            results[s][i] = Array(t, n)
        end
    end

    scenario = 0

    # Initialise the objective storage
    obj = 0.

    # For all stages
    for stage=1:m.stages
        # realise on scenario
        sp = m.stage_problems[stage, markov[stage]]
        scenario = load_scenario!(m, sp)
        data = stagedata(sp)
        for (key, series) in kwargs
            @assert haskey(data.scenario_constraint_names, key)
            cidx = data.scenario_constraint_names[key]
            JuMP.setRHS(data.scenario_constraints[cidx][1], series[stage])
        end

        # solve
        solve!(sp)

        # Add objective (stage profit only)
        obj += get_true_value(sp)

        # Save results if necesary
        store_results!(results, vars, sp, stage, 1, markov[stage], scenario)
        # pass forward if necessary
        if stage < m.stages
            pass_states_forward!(m, stage, markov[stage])
        end
    end

    results
end
