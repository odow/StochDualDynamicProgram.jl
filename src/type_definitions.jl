#  Copyright 2017, Oscar Dowson

typealias LinearConstraint JuMP.ConstraintRef{JuMP.Model, JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64, JuMP.Variable}}}

abstract Sense
immutable Minimisation <: Sense end
immutable Maximisation <: Sense end

"""
    SDDP State Variable
"""
immutable StateVariable
    x::JuMP.Variable
    con::LinearConstraint
end

"""
    Scenario constraint
"""
immutable Scenario{N}
    con::LinearConstraint
    values::Vector{Float64}
end

"""
    All scenarios
"""
type Scenarios{N}
    con::Vector{Scenario{N}}
    prob::Vector{Float64}
end
Scenarios(probabilities) = Scenarios{length(probabilities)}(Scenario{length(probabilities)}[], probabilities)

immutable DiscreteDistribution{T}
    values::Vector{T}
    support::Vector{Float64}
end
DiscreteDistribution(x) = DiscreteDistribution(x, ones(length(x)) / length(x))

"""
    Price scenario information
"""
type PriceScenarios{T}
    ribs::Vector{Float64}
    noises::DiscreteDistribution{T}
    dynamics::Function
    objective::Function
    lastprice::Float64
end

"""
    Abstract type for dispatching the cut function
"""
abstract AbstractRiskMeasure

"""
    Normal old expectation
"""
immutable Expectation end

"""
    Nested CV@R
        λE[x] + (1-λ)CV@R(1-α)(x)
"""
immutable NestedCVaR
    alpha::Float64
    lambda::Float64
end

setriskmeasure!(m::JuMP.Model, riskmeasure::AbstractRiskMeasure) = (ext(m).riskmeasure = riskmeasure)

"""
    The extension type for the JuMP subproblems
"""
type SubproblemExt{S, P, RM<:AbstractRiskMeasure}
    states::Vector{StateVar}
    scenarios::Scenarios{S}
    pricescenarios::PriceScenarios{P}
    riskmeasure::RM
    stageobjective::JuMP.GenericAffExpr{Float64, JuMP.Variable}
end
function addsubproblemdata!(m, scenario_probabilities, pricescenarios::Int, riskmeasure)
    m.ext[:subproblem] = SubproblemExt(
        StateVar[],
        Scenarios(scenario_probabilities),

        riskmeasure,
        JuMP.AffExpr()
    )
end

ext(m::JuMP.Model) = m.ext[:subproblem]::SubproblemExt

"""
    Θ (</>)= intercept + coefficients' x
"""
immutable Cut{N}
    intercept::Float64
    coefficients::Tuple{Vararg{Float64, N}}
end

type CutStorage{N}
    cuts::Vector{Cut{N}}
    states_dominant::Vector{Int} # states_dominant[i] = number of states in statesvisited that cut i is dominant
    statesvisited::Vector{Tuple{Vararg{Float64, N}}} # vector of state tuples visited
    best_cut_index::Vector{Int} # best_cut_index[i] = index of cut in cuts that is the dominant cut at states_dominant[i]
    best_bound::Vector{Float64} # best_bound[i] = best objective bound at statesvisited[i]
end
CutStorage(N) = CutStorage{N}(Cut{N}[], Int[], Tuple{Vararg{Float64, N}}[], Int[], Float64[])

"""
    The main type that holds the Stochastic Dual Dynamic Programming model
"""
type SDDPModel{N}
    # subproblems
    stageproblems::Vector{Vector{JuMP.Model}}
    # corresponding cut storage
    cutstorage::Vector{Vector{CutStorage{N}}}

    # markov transition matrices
    transition::Vector{Array{Float64, 2}}

    buildsubproblem!::Function
end

"""
    Some accessor functions
"""
stageproblem(m::SDDPModel, t::Int, i::Int) = m.stageproblems[t][i]
stageproblem(m::SDDPModel, t::Int) = stageproblem(m, t, 1)
cutstorage(m::SDDPModel, t::Int, i::Int) = m.cutstorage[t][i]
cutstorage(m::SDDPModel, t::Int) = cutstorage(m, t, 1)

function SDDPModel(buildsubproblem!::Function;
    sense::Symbol=:Minimisation,
    stages::Int = 1,
    transition = [1]',
    scenarios = [1.],
    pricescenarios = 1
    )

    stageproblems = Vector{JuMP.Model}[]
    for t=1:stages
        push!(stageproblems, JuMP.Model[])
        for i=1:nummarkovstates(transition, t)
            m = JuMP.Model()
            addsubproblemdata!(m,
                getvec(scenarios, t, i),
                getel(pricescenarios, t, i)
            )
            buildsubproblem!(m, t, i)
            push!(stageproblems[t], m)
        end
    end

    N = length(ext(stageproblems[1][1]).state_vars)
    SDDPModel(
        stageproblems,
        [CutStorate{N}[CutStorage(N) for i=1:nummarkovstates(transition, t)] for t=1:stages]
        transition,
        buildsubproblem!
    )

end

nummarkovstates(T::Array{Float64, 2}, t::Int) = size(T, 1)
nummarkovstates(T::Vector{Array{Float64, 2}}, t::Int) = size(T[t], 1)

getvec{T<:AbstractVector}(x::T, t::Int, i::Int) = x
getvec{T<:AbstractVector}(x::Vector{T}, t::Int, i::Int) = x[t]
getvec{T<:AbstractVector}(x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
getel{T}(x::T, t::Int, i::Int) = x
getel{T}(x::Vector{T}, t::Int, i::Int) = x[t]
getel{T}(x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
