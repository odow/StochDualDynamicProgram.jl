#  Copyright 2017, Oscar Dowson

function stageobjective!(m::JuMP.Model, stageobj::JuMP.GenericAffExpr)
    @assert length(ext(m).theta) == 0
    theta = @variable(m)
    push!(getribs(m), Rib(0.0, theta))
    _setobjective!(m, theta + stageobj)
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
            getribs(sp),
            Rib(r, @variable(sp))
        )
    end
    for i in 1:length(noises)
        push!(
            getpricescenarios(sp),
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
