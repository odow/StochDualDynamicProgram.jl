#  Copyright 2017, Oscar Dowson

"""
    Normal old expectation
"""
immutable Expectation <: AbstractRiskMeasure end

function cutgenerator(ex::Expectation, m, x, pi, theta, prob)
    @assert length(pi) == length(theta) == length(prob)
    intercept = (theta[1] - dot(pi[1], x))*prob[1]
    coefficients = pi[1] * prob[1]
    @inbounds for i=2:length(prob)
        intercept += (theta[i] - dot(pi[i], x))*prob[i]
        coefficients += pi[i] * prob[i]
    end
    return Cut(intercept, coefficients)
end
