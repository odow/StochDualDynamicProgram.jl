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
JuMP.getdual(x::StateVariable) = getdual(x.con)

"""
    Scenario constraint
"""
immutable Scenario
    probability::Float64
    con::Vector{LinearConstraint}
    values::Vector{Float64}
end
Scenario(prob::Float64) = Scenario(prob, LinearConstraint[], Float64[])

immutable DiscreteDistribution{T}
    values::Vector{T}
    support::Vector{Float64}
end
DiscreteDistribution(x) = DiscreteDistribution(x, ones(length(x)) / length(x))

immutable Rib
    price::Float64
    theta::JuMP.Variable
end

"""
    Price scenario information
"""

type PriceScenarios
    noises::Vector
    probability::Vector{Float64}
    dynamics::Function
    objective::Function
end
PriceScenarios() = PriceScenarios(Any[], Float64[], ()->(), ()->())

"""
    Abstract type for dispatching the cut function
"""
abstract AbstractRiskMeasure

"""
    Normal old expectation
"""
immutable Expectation <: AbstractRiskMeasure end

"""
    Nested CV@R
        λE[x] + (1-λ)CV@R(1-α)(x)
"""
immutable NestedCVaR <: AbstractRiskMeasure
    alpha::Float64
    lambda::Float64
    storage::Vector{Float64}
end
NestedCVaR(alpha, lambda) = NestedCVaR(alpha, lambda, Float64[])

"""
    The extension type for the JuMP subproblems
"""
type SubproblemExt{T}
    states::Vector{StateVar}
    scenarios::Vector{Scenario}
    pricescenarios::PriceScenarios{T}
    theta::Vector{Rib}
end
function Subproblem()
    m = JuMP.Model()
    m.ext[:subproblem] = SubproblemExt(
        StateVar[],
        Scenario[],
        PriceScenarios(),
        JuMP.AffExpr(),
        Rib[]
    )
    m
end

ext(m::JuMP.Model) = m.ext[:subproblem]::SubproblemExt

"""
    Θ (</>)= intercept + coefficients' x
"""
immutable Cut
    intercept::Float64
    coefficients::Vector{Float64}
end
Cut(intercept, coefficients::Vector) = Cut(intercept, coefficients)

# There is a cut oracle for every subproblem (and every price state).
abstract CutOracle
typealias OracleStore{T} Vector{Vector{Vector{T}}}
function OracleStore(T::CutOracle, sense, problem_size)
    [
        [
            [
                initialise(T, sense stage, markovstate, pricestate)
            for pricestate in 1:num_price_states]
        for (markovstate, num_price_states) in enumerate(markov_states)]
    for (stage, markov_states) in enumerate(problem_size)]
end

"""
    The main type that holds the Stochastic Dual Dynamic Programming model
"""
type SDDPModel{C, R}
    # subproblems
    stageproblems::Vector{Vector{JuMP.Model}}
    # Cut oracle
    cutoracle::OracleStore{C}
    # list of states visited in each stage
    statesvisited::Vector{Vector{Vector{Float64}}}
    # Risk Measure
    riskmeasure::R

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

numstages(m::SDDPModel) = length(m.stageproblems)
nummarkovstates(m::SDDPModel, stage::Int) = length(m.stageproblems[stage])
numpricescenarios(m::JuMP.Model) = length(ext(m).pricescenarios.ribs)
numpricescenarios(m::SDDPModel, stage::Int, markovstate::Int) = numpricescenarios(stageproblem(stage, markovstate))
numpricescenarios(m::SDDPModel, stage::Int) = numpricescenarios(m, stage, 1)

numscenarios(m::JuMP.Model) = length(ext(m).scenarios.prob)
numscenarios(m::SDDPModel, stage::Int, markovstate::Int) = numscenarios(stageproblem(stage, markovstate))
numscenarios(m::SDDPModel, stage::Int) = numscenarios(m, stage, 1)

function initialisescenaros!(m::JuMP.Model, scenario_probabilities::Vector{Float64})
    if (sum(scenarios) - 1.0) > 1e-6
        error("Sum of scenario probabilities must sum to 1. You have given the vector $(scenario_probabilities) which sums to $(sum(scenario_probabilities))")
    end
    for p in scenario_probabilities
        push!(ext(m).scenarios, Scenario(p))
    end
end

function SDDPModel(buildsubproblem!::Function;
    sense::Symbol=:Min,
    stages::Int = 1,
    transition = [1]',
    riskmeasure = Expectation(),
    cutoracle   = DefaultCutOracle(),
    scenarios = 0
    )

    stageproblems = Vector{JuMP.Model}[]
    problem_size = Vector{Int}[]
    for t=1:stages
        push!(problem_size, Int[])
        push!(stageproblems, JuMP.Model[])
        for i=1:nummarkovstates(transition, t)
            m = Subproblem()
            JuMP.setobjectivesense(m, sense)
            buildsubproblem!(m, t, i)
            initialisescenarios!(m, getscenariovec(scenarios, t, i))
            #
            push!(stageproblems[t], m)
            push!(problem_size[t], numpricescenarios(m, t, i))
        end
    end
    SDDPModel(
        stageproblems,
        OracleStore(cutoracle, getsense(sense), problem_size),
        riskmeasure,
        transition,
        buildsubproblem!
    )

end

function getsense(x::Symbol)
    if x == :Min
        return Minimisation
    elseif x == :Max
        return Maximisation
    else
        error("The sense $(x) must be one of [:Min, :Max]")
    end
end

nummarkovstates(T::Array{Float64, 2}, t::Int) = size(T, 1)
nummarkovstates(T::Vector{Array{Float64, 2}}, t::Int) = size(T[t], 1)
getscenariovec(x::Int, t::Int, i::Int) = ones(x) / x
getscenariovec{T<:AbstractVector}(x::T, t::Int, i::Int) = x
getscenariovec{T<:AbstractVector}(x::Vector{T}, t::Int, i::Int) = x[t]
getscenariovec{T<:AbstractVector}(x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
getel{T}(x::T, t::Int, i::Int) = x
getel{T}(x::Vector{T}, t::Int, i::Int) = x[t]
getel{T}(x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
