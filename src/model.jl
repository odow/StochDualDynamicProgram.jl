function initialisescenarios!(m::JuMP.Model, scenario_probabilities::Vector{Float64})
    if (sum(scenario_probabilities) - 1.0) > 1e-6
        error("Sum of scenario probabilities must sum to 1. You have given the vector $(scenario_probabilities) which sums to $(sum(scenario_probabilities))")
    end
    for p in scenario_probabilities
        push!(ext(m).scenarios, Scenario(p))
    end
end

function SDDPModel(buildsubproblem!::Function;
    sense::Symbol=:Min,
    stages::Int = 1,
    transition = [1.0]',
    riskmeasure = Expectation(),
    cutoracle   = DefaultCutOracle(),
    scenarios = 0
    )

    stageproblems = Vector{JuMP.Model}[]
    statesvisited = Vector{Vector{Float64}}[]
    problem_size = Vector{Int}[]

    tmp_storage_size = 0
    max_size = 0
    for t=1:stages

        push!(problem_size, Int[])
        push!(stageproblems, JuMP.Model[])
        push!(statesvisited, Vector{Float64}[])
        max_size = 0
        for i=1:nummarkovstates(transition, t)
            m = Subproblem()
            JuMP.setobjectivesense(m, sense)
            initialisescenarios!(m, getscenariovec(scenarios, t, i))
            buildsubproblem!(m, t, i)

            push!(stageproblems[t], m)
            push!(problem_size[t], numpriceribs(m))

            max_size += numscenarios(m) * numpricescenarios(m)
        end
        if max_size > tmp_storage_size
            tmp_storage_size = max_size
        end
    end
    tmp_storage = TmpStorage(tmp_storage_size, numstates(stageproblems[1][1]))
    SDDPModel{getsense(sense), typeof(cutoracle), typeof(riskmeasure), typeof(transition)}(
        stageproblems,
        OracleStore(cutoracle, getsense(sense), problem_size),
        statesvisited,
        riskmeasure,
        transition,
        buildsubproblem!,
        tmp_storage
    )

end
