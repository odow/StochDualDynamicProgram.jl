function initialisescenarios!(m::JuMP.Model, scenario_probabilities::Vector{Float64})
    if (sum(scenario_probabilities) - 1.0) > 1e-6
        error("Sum of scenario probabilities must sum to 1. You have given the vector $(scenario_probabilities) which sums to $(sum(scenario_probabilities))")
    end
    for p in scenario_probabilities
        push!(ext(m).scenarios, Scenario(p))
    end
end

validatetransition(transition::Array{Float64, 2}, initialmarkovstate) = nothing
function validatetransition(transition::Vector{Array{Float64, 2}}, stages::Int)
    if length(transition) != (stages - 1)
            error("Incorrect markov transition supplied. You have chosen to give a vector of transition matrices with a specified markov state for the first stage, so you should have a vector that contains (stages - 1) = $(stages-1) markov transition matrices. However the current vector has $(length(transition)).")
    end
end
initialmarkovprobability(x, n) = ones(n) / n
function initialmarkovprobability(x::Vector, n)
    @assert length(x) == n
    @assert (sum(x) - 1.0) < 1e-6
    return x
end
function SDDPModel(buildsubproblem!::Function;
    sense::Symbol=:Min,
    stages::Int = 1,
    transition = [1.0]',
    riskmeasure = Expectation(),
    cutoracle   = DefaultCutOracle(),
    scenarios = 0,
    initialmarkovstate = nothing,
    firstprice = NaN
    )
    validatetransition(transition, stages)

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
    SDDPModel{getsense(sense), typeof(cutoracle), typeof(riskmeasure), typeof(transition)}(
        stageproblems,
        OracleStore(cutoracle, getsense(sense), problem_size),
        statesvisited,
        riskmeasure,
        transition,
        initialmarkovprobability(initialmarkovstate, nummarkovstates(transition, 1)),
        firstprice,
        buildsubproblem!,
        [BackwardStorage(numstates(stageproblems[1][1])) for i=1:tmp_storage_size],
        ForwardStorage[]
    )

end
