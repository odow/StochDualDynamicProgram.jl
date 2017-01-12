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

    m = SDDPModel(buildsubproblem!)

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

function priceprocess!(sp::JuMP.Model,
    ribs::AbstractVector{Float64},
    dynamics::Function,
    noises::AbstractVector,
    probability::AbstractVector{Float64},
    objective::Function
    )
    @assert length(noises) == length(probability)
    for r in ribs
        push!(
            ribs(sp),
            Rib(r, @variable(sp))
        )
    end
    for i in 1:length(noises)
        push!(
            pricescenarios(sp),
            PriceScenario(probability[i], createpricefunction!(dynamics, objective, noises[i]))
        )
    end

end
priceprocess!(sp::JuMP.Model,
    ribs::AbstractVector{Float64},
    dynamics::Function,
    noises::AbstractVector,
    # probability::AbstractVector{Float64},
    objective::Function
    ) = priceprocess!(sp, ribs, dynamics, noises, ones(length(noises)) / length(noises), objective)

function createpricefunction!(dynamics::Function, objective::Function, noise)
    (p) -> (
        newp = dynamics(p, noise);
        (newp, objective(newp))
    )
end
