#  Copyright 2017, Oscar Dowson

typealias LinearConstraint JuMP.ConstraintRef{JuMP.Model, JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64, JuMP.Variable}}}

abstract AbstractSense
immutable Minimisation <: AbstractSense end
immutable Maximisation <: AbstractSense end

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
    con::LinearConstraint
    values::Vector{Float64}
end

"""
    All scenarios
"""
type Scenarios
    con::Vector{Scenario}
    prob::Vector{Float64}
end
Scenarios() = Scenarios(Scenario[], Float64[])

immutable DiscreteDistribution{T}
    values::Vector{T}
    support::Vector{Float64}
end
DiscreteDistribution(x) = DiscreteDistribution(x, ones(length(x)) / length(x))

"""
    Price scenario information
"""
type PriceScenarios
    ribs::Vector{Float64}
    noises::DiscreteDistribution # TODO: not type stable
    dynamics::Function
    objective::Function
    lastprice::Float64
end
PriceScenarios() = PriceScenarios(Float64[] DiscreteDistribution([0.], [1.]), ()->(), ()->(), 0.)

"""
    Abstract type for dispatching the cut function
"""
abstract CutOracle

"""
    Normal old expectation
"""
immutable Expectation <: CutOracle end

"""
    Nested CV@R
        λE[x] + (1-λ)CV@R(1-α)(x)
"""
immutable NestedCVaR <: CutOracle
    alpha::Float64
    lambda::Float64
    storage::Vector{Float64}
end
NestedCVaR(alpha, lambda) = NestedCVaR(alpha, lambda, Float64[])

# setriskmeasure!(m::JuMP.Model, riskmeasure::AbstractRiskMeasure) = (ext(m).riskmeasure = riskmeasure)

"""
    The extension type for the JuMP subproblems
"""
type SubproblemExt{RM<:AbstractRiskMeasure}
    states::Vector{StateVar}
    scenarios::Scenarios
    pricescenarios::PriceScenarios
    riskmeasure::RM
    stageobjective::JuMP.GenericAffExpr{Float64, JuMP.Variable}
    theta::JuMP.Variable
end
function addsubproblemdata!(m, scenario_probabilities, riskmeasure)
    m.ext[:subproblem] = SubproblemExt(
        StateVar[],
        Scenarios(),
        PriceScenarios(),
        riskmeasure,
        JuMP.AffExpr(),
        @variable(m)
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

type CutStorageInner{N}
    cuts::Vector{Cut{N}}
    states_dominant::Vector{Int} # states_dominant[i] = number of states in statesvisited that cut i is dominant
    best_cut_index::Vector{Int} # best_cut_index[i] = index of cut in cuts that is the dominant cut at states_dominant[i]
    best_bound::Vector{Float64} # best_bound[i] = best objective bound at statesvisited[i]
end
CutStorageInner(N) = CutStorageInner{N}(Cut{N}[], Int[], Int[], Float64[])

type CutStorage{N}
    cutstorage::CutStorageInner{N}
    statesvisited::Vector{Tuple{Vararg{Float64, N}}} # vector of state tuples visited
end
CutStorage(N) = CutStorage(CutStorageInner(N), Tuple{Vararg{Float64, N}}[])

"""
    The main type that holds the Stochastic Dual Dynamic Programming model
"""
type SDDPModel{N}
    # subproblems
    stageproblems::Vector{Vector{JuMP.Model}}
    # corresponding cut storage cutstoreage[stage][markovstate][rib]
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
    sense::Symbol=:Min,
    stages::Int = 1,
    transition = [1]',
    riskmeasure = Expectation()
    )

    stageproblems = Vector{JuMP.Model}[]
    for t=1:stages
        push!(stageproblems, JuMP.Model[])
        for i=1:nummarkovstates(transition, t)
            m = JuMP.Model()
            addsubproblemdata!(m,
                getel(riskmeasure, t, i)
            )
            buildsubproblem!(m, t, i)
            JuMP.setobjective(m, sense, ext(m).theta + ext(m).stageobjective)
            push!(stageproblems[t], m)
        end
    end

    N = length(ext(stageproblems[1][1]).state_vars)
    SDDPModel(
        stageproblems,
        [CutStorage{N}[CutStorage(N) for i=1:nummarkovstates(transition, t)] for t=1:stages]
        transition,
        buildsubproblem!
    )

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
