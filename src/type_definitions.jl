#  Copyright 2017, Oscar Dowson

typealias LinearConstraint JuMP.ConstraintRef{JuMP.Model, JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64, JuMP.Variable}}}

typealias Sense Union{Type{Val{:Min}}, Type{Val{:Max}}}
typealias Minimisation Val{:Min}
typealias Maximisation Val{:Max}

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
immutable Scenario
    probability::Float64
    con::Vector{LinearConstraint}
    values::Vector{Float64}
end
Scenario(prob::Float64) = Scenario(prob, LinearConstraint[], Float64[])

typealias DiscreteDistribution{T} Vector{Tuple{T, Float64}}
DiscreteDistribution{T<:Real}(x::AbstractVector{T}) = map(i->(i, 1.0/length(x)), x)
function DiscreteDistribution{T, T2<:Real}(x::AbstractVector{T}, y::AbstractVector{T2})
    @assert length(x) == length(y)
    map((i, j) -> (i, j), x, y)
end
immutable Rib
    price::Float64
    theta::JuMP.Variable
end

"""
    Price scenario information
"""

type PriceScenario{T}
    noise::T
    probability::Float64
    dynamics::Function
    objective::Function
end
PriceScenario() = PriceScenario(0.0, 1.0, ()->(), ()->())

"""
    Abstract type for dispatching the cut function
"""
abstract AbstractRiskMeasure

"""
    The extension type for the JuMP subproblems
"""
type SubproblemExt
    states::Vector{StateVariable}
    scenarios::Vector{Scenario}
    pricescenarios::Vector{PriceScenario}
    theta::Vector{Rib}
end
function Subproblem()
    m = JuMP.Model()
    m.ext[:subproblem] = SubproblemExt(
        StateVariable[],
        Scenario[],
        PriceScenario[],
        Rib[]
    )
    m
end

"""
    Θ (</>)= intercept + coefficients' x
"""
immutable Cut
    intercept::Float64
    coefficients::Vector{Float64}
end
Cut(intercept, coefficients::Vector) = Cut(intercept, coefficients)

immutable CutContainer
    cut::Cut
    stage::Int
    markov::Int
    rib::Int
end

# There is a cut oracle for every subproblem (and every price state).
abstract CutOracle
typealias OracleStore{T} Vector{Vector{Vector{T}}}
function OracleStore(T::CutOracle, sense, problem_size)
    [
        [
            [
                initialise(T, sense, stage, markovstate, pricestate)
            for pricestate in 1:num_price_states]
        for (markovstate, num_price_states) in enumerate(markov_states)]
    for (stage, markov_states) in enumerate(problem_size)]
end

type BackwardStorage
    obj::Float64
    probability::Float64
    pi::Vector{Float64}
end

BackwardStorage(xn) = BackwardStorage(0.0, 0.0, zeros(Float64, xn))

type ForwardStorage
    markov::Int
    price::Float64
    state::Vector{Float64}
end


function ForwardStorage(stages, states)
    [ForwardStorage(0, 0.0, zeros(Float64, states)) for i=1:stages]
end

"""
    The main type that holds the Stochastic Dual Dynamic Programming model
"""
type SDDPModel{S, C, R, F1, F2, T}
    # subproblems
    stageproblems::Vector{Vector{JuMP.Model}}
    # Cut oracle
    cutoracle::OracleStore{C}
    # list of states visited in each stage
    statesvisited::Vector{Vector{Vector{Float64}}}
    # Risk Measure
    riskmeasure::R

    forwardsampler::F1 # forward sampler

    # markov transition matrices
    transition::T
    initialmarkovprobability::Vector{Float64}

    initialprice::Float64

    buildsubproblem!::Function

    backwardstorage::Vector{BackwardStorage}
    forwardstorage::Vector{Vector{ForwardStorage}}
end
