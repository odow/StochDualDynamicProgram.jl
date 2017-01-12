#  Copyright 2017, Oscar Dowson

"""
    This function assembles a new cut using the following inputs
    + measure::AbstractRiskMeasure - used to dispatch
    + m::SDDPModel
    + x::Vector{Vector{Float64}} - a vector of vector of state values for each scenario
    + pi::Vector{Vector{Float64}} - a vector of vector of dual values for each scenario
    + theta::Vector{Float64} - a vector of objective values for each scenario
    + prob::Vector{Float64} - the probability support of the scenarios. Should sum to one
    + stage::Int - the index of the stage
    + markov::Int - the index of the markov state
"""
cutgenerator(measure::AbstractRiskMeasure, m::SDDPModel, x, pi, theta, prob, stage, markov) =
    cutgenerator(measure::AbstractRiskMeasure, m::SDDPModel, x, pi, theta, prob)

cutgenerator(measure::AbstractRiskMeasure, m::SDDPModel, x, pi, theta, prob) = error("""
    You need to overload a `cutgenerator` method for the measure of type $(typeof(measure)).
    This could be the method including the stage and markov index
        cutgenerator(measure::AbstractRiskMeasure, sense::Sense, x, pi, theta, prob, stage, markov)
    or
        cutgenerator(measure::AbstractRiskMeasure, sense::Sense, x, pi, theta, prob)
""")

include("riskmeasures/expectation.jl")
include("riskmeasures/nestedcvar.jl")
