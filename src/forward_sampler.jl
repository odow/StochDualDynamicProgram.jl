# Copyright 2017, Oscar Dowson

# This file here is really only extra for experts (and because @lkapelevich asked for it)
function sample(x::Vector, w=ones(length(x))/length(x))
    @assert length(w) == length(x)
    @assert (sum(x) - 1.0) <= 1e-6
    r = rand(Float64)
    @inbounds for i in 1:length(x)
        r -= w[i]
        if r < 1e-6
            return x[i]
        end
    end
end

abstract AbstractForwardSampler

samplemarkov(sampler::AbstractForwardSampler, transitionmatrix::Array{Float64, 2}, stage, lastmarkov) = error("""
    You must define a samplemarkov method for your forward sampler $(typeof(sampler)).
""")

samplescenario(sampler::AbstractForwardSampler, scenarios::Vector{Scenario}, stage, lastmarkov) = error("""
    You must define a samplescenario method for your forward sampler $(typeof(sampler)).
""")

samplepricescenario(sampler::AbstractForwardSampler, pricescenarios::Vector{PriceScenario}, stage, lastmarkov) = error("""
    You must define a samplepricescenario method for your forward sampler $(typeof(sampler)).
""")


"""
    Sample the forward pass with the a uniform probability across all outcomes.
"""
immutable UniformSampler <: ForwardSampler end

samplemarkov(::UniformSampler, transitionmatrix::Array{Float64, 2}, stage, lastmarkov) = sample(1:size(transitionmatrix, 2))
samplescenario(::UniformSampler, scenarios::Vector{Scenario}, stage, markov) = sample(scenarios)
samplepricescenario(::UniformSampler, pricescenarios::Vector{PriceScenario}, stage, markov) = sample(pricescenarios)

"""
    Sample the forward pass with the same probabilities as those used to simulate
    the policy.
"""
immutable DefaultSampler <: ForwardSampler end

samplemarkov(::UniformSampler, transitionmatrix::Array{Float64, 2}, stage, lastmarkov) = sample(1:size(transitionmatrix, 2), transitionmatrix[lastmarkov, :])

function scenariosampler{T<:Union{Scenario, PriceScenario}}(::DefaultSampler, scenarios::Vector{T}, stage, markov)
    r = rand()
    for scenario in scenarios
        r -= probability(scenario)
        if r <= 0.0
            return scenario
        end
    end
end

samplescenario(sampler::DefaultSampler, scenarios::Vector{Scenario}, stage, markov) = scenariosampler(sampler, scenarios, stage, markov)

samplepricescenario(::DefaultSampler, pricescenarios::Vector{PriceScenario}, stage, markov) = scenariosampler(sampler, pricescenarios, stage, markov)
