# Copyright 2017, Oscar Dowson

function ExtendedJuMPModel()
    sp = JuMP.Model()
    sp.ext[:subproblem] = SubproblemExtension()
    sp
end
ext(m::JuMP.Model) = m.ext[:subproblem]::SubproblemExtension


function SDDPModel(buildsubproblem!::Function;
    stages::Int = 1,
    sense::Symbol = :Min,
    markovstates = 1,
    transition = nothing,
    scenarios = nothing,
    kwargs...)

    m = SDDPModel(getsense(sense), buildsubproblem!)

    for t in 1:stages
        stage = Stage()
        for i in 1:numberofmarkovstates(markovstates, t)
            sp = ExtendedJuMPModel()
            JuMP.setobjectivesense(sp, sense)
            settransition!(sp, transition, markovstates, stages, t, i)
            setscenarios!(sp, scenarios, t, i)
            buildsubproblem!(sp, t, i)
            push!(stage.arr, Subproblem(sp, DefaultCutOracle(), Expectation()))
        end
        push!(m.stageproblems, stage)
    end
    m
end
