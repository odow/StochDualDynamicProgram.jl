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

"""
    Price scenario information
"""
type PriceScenarios{T<:AbstractVector}
    ribs::Vector{Float64}
    noises::T
    probabilitysupport::Vector{Float64}
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

ext(m::JuMP.Model) = m.ext::SubproblemExt

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
    solver::MathProgBase.AbstractMathProgSolver
end

"""
    Some accessor functions
"""
stageproblem(m::SDDPModel, t::Int, i::Int) = m.stageproblems[t][i]
stageproblem(m::SDDPModel, t::Int) = stageproblem(m, t, 1)
cutstorage(m::SDDPModel, t::Int, i::Int) = m.cutstorage[t][i]
cutstorage(m::SDDPModel, t::Int) = cutstorage(m, t, 1)
