#  Copyright 2017, Oscar Dowson

"""
    This function assembles a new cut using the following inputs
    + measure::AbstractRiskMeasure - used to dispatch
    + sense::Sense - either Maximum or Minimum
    + x::Vector{Vector{Float64}} - a vector of vector of state values for each scenario
    + pi::Vector{Vector{Float64}} - a vector of vector of dual values for each scenario
    + theta::Vector{Float64} - a vector of objective values for each scenario
    + prob::Vector{Float64} - the probability support of the scenarios. Should sum to one
    + stage::Int - the index of the stage
    + markov::Int - the index of the markov state
"""
cutgenerator(measure::AbstractRiskMeasure, sense::Sense, x, pi, theta, prob, stage, markov) =
    cutgenerator(measure::AbstractRiskMeasure, sense::Sense, x, pi, theta, prob)

cutgenerator(measure::AbstractRiskMeasure, sense::Sense, x, pi, theta, prob) = error("""
    You need to overload a `cutgenerator` method for the measure of type $(typeof(measure)).
    This could be the method including the stage and markov index
        cutgenerator(measure::AbstractRiskMeasure, sense::Sense, x, pi, theta, prob, stage, markov)
    or
        cutgenerator(measure::AbstractRiskMeasure, sense::Sense, x, pi, theta, prob)
""")


"""
    The expectation risk measure
"""
function cutgenerator(ex::Expectation, sense::Sense, x, pi, theta, prob)
    @assert length(pi) == length(theta) == length(prob)
    intercept = theta[1] - dot(pi[1], x)
    coefficients = pi[1] * prob[1]
    @inbounds for i=2:length(prob)
        intercept += theta[i] - dot(pi[i], x)
        coefficients += pi[i] * prob[i]
    end
    return Cut(intercept, coefficients)
end

# λ * E[x] + (1 - λ) * CVaR(β)[x]
const expectation = Expectation()
_sortperm(::Type{Maximisation}, x) = sortperm(x, rev=false)
_sortperm(::Type{Minimisation}, x) = sortperm(x, rev=true)
function calculateCVaRprobabilities!(newprob, sense, oldprob, theta, lambda, beta::Float64)
    @assert length(newprob) >= length(oldprob)
    quantile_collected = 0.0
    cvarprob = 0.0
    cache = (1 - lambda) / beta                # cache this to save some operations
    @inbounds for i in _sortperm(sense, theta) # For each scenario in order
        newprob[i] = lambda * oldprob[i]       # expectation contribution
        if quantile_collected <  beta          # We haven't collected the beta quantile
            cvarprob = min(oldprob[i], beta-quantile_collected)
            newprob[i] += cache * cvarprob     # risk averse contribution
            quantile_collected += cvarprob     # Update total quantile collected
        end
    end
end

"""
    A weighted combination of Expectation and CVaR
"""
function cutgenerator(cvar::NestedCVaR, sense::Sense, x, pi, theta, prob)
    @assert length(pi) == length(theta) == length(prob)
    if length(cvar.storage) < length(prob)
        append!(cvar.storage, zeros(length(prob) - length(cvar.storage)))
    end
    calculateCVaRprobabilities!(cvar.storage, sense, prob, x, cvar.lambda, cvar.beta)

    cutgenerator(expectation, x, pi, theta, cvar.storage[1:length(prob)])
end
