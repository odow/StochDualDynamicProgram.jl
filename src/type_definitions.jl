#  Copyright 2016, Oscar Dowson

"""
    The main type that holds the Stochastic Dual Dynamic Programming model
"""
type SDDPModel
    stageproblems::Vector{Vector{JuMP.Model}}
end

"""
    Some accessor functions
"""
stageproblem(m::SDDPModel, t::Int, i::Int) = m.stageproblems[t][i]
stageproblem(m::SDDPModel, t::Int) = stageproblem(m, t, 1)

"""
    SDDP State Variable
"""
immutable StateVariable{C}
    x::JuMP.Variable
    con::C
end

"""
    Scenario constraint
"""
immutable Scenario{C, N}
    con::C
    values::Vector{Float64}
end

"""
    All scenarios
"""
type Scenarios{C, N}
    con::Vector{Scenario{C, N}}
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
type SubproblemExt{C1, C2, N, T, RM<:AbstractRiskMeasure}
    states::Vector{StateVar{C1}}
    scenarios::Scenarios{C2, N}
    pricescenarios::PriceScenarios{T}
    riskmeasure::RM
    stageobjective::JuMP.GenericAffExpr{Float64, JuMP.Variable}
end

ext(m::JuMP.Model) = m.ext::SubproblemExt
