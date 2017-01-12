#  Copyright 2017, Oscar Dowson

"""
    Nested CV@R
        λE[x] + (1-λ)CV@R(1-α)(x)
"""
immutable NestedCVaR <: AbstractRiskMeasure
    beta::Float64
    lambda::Float64
    storage::Vector{Float64}
end

function NestedCVaR(beta, lambda)
    checkzerotoone(beta)
    checkzerotoone(lambda)
    NestedCVaR(beta, lambda, Float64[])
end
NestedCVaR(;beta=1, lambda=1) = NestedCVaR(beta, lambda)

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

function cutgenerator(cvar::NestedCVaR, m, x, pi, theta, prob)
    @assert length(pi) == length(theta) == length(prob)
    if length(cvar.storage) < length(prob)
        append!(cvar.storage, zeros(length(prob) - length(cvar.storage)))
    end
    calculateCVaRprobabilities!(cvar.storage, sense, prob, theta, cvar.lambda, cvar.beta)
    cutgenerator(expectation, m, x, pi, theta, view(cvar.storage, 1:length(prob)))
end
