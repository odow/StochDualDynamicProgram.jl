#  Copyright 2017, Oscar Dowson

"""
    @stageobjective(sp, ex)

Define the stage profit for subproblem `sp`.

Arguments:

    sp    the subproblem
    ex    a JuMP expression for the costs accrued in the subproblem

Usage:

    @stageobjective(sp, x + y)
"""
macro stageobjective(m, ex)
    m = esc(m)
    quote
        depwarn("""
            The macro @stageobjective has been deprecated. Use the function stageobjective! from now on. It has the same syntax as before. i.e.:

                @stageobjective(m, 2x + y)

                becomes

                stageobjective!(m, 2x + y)
                """)
        stageobjective!($m, $(esc(ex)))
    end
end

function stageobjective!(m::JuMP.Model, stageobj::JuMP.GenericAffExpr)
    @assert length(ext(m).theta) == 0
    theta = JuMP.@variable(m)
    push!(ribs(m), Rib(0.0, theta))
    _setobjective!(m, theta + stageobj)
end

"""
    objectivescenario!(sp,                 # subproblem
        rib_locations = 0:10,              # outgoing price
        noises        = rand(10),          # noises
        dynamics      = (p, w) -> (p + w), # dynamics
        objective     = (p) -> (p * x)     # Objective can use p0 as a parameter, x as a variable
    )
"""
objectivescenario!{T<:Real}(sp::JuMP.Model;
    rib_locations::AbstractVector{T}=Float64[],
    noises=Float64[],
    dynamics::Function=()->(),
    objective::Function=()->()
    ) = objectivescenario!(sp,
        rib_locations,
        noises,
        dynamics,
        objective
    )
objectivescenario!(sp::JuMP.Model, discretisation::AbstractVector, noises::AbstractVector, dynamics::Function, objective::Function) = objectivescenario!(sp, discretisation, DiscreteDistribution(noises), dynamics, objective)

function objectivescenario!(sp::JuMP.Model, discretisation::AbstractVector, noises::DiscreteDistribution, dynamics::Function, objective::Function)
    ex = ext(sp)
    if length(ex.theta) > 0
        error("""
        You can only define one objectivescenario for each subproblem.

        Do you have a call to stageobjective! as well?
            """)
    end
    for rib in discretisation
        push!(ex.theta, Rib(rib, @variable(sp)))
    end
    for (val, prob) in noises
        push!(ex.pricescenarios, PriceScenario(val, prob, dynamics, objective))
    end
end
