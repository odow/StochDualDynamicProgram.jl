# Copyright 2017, Oscar Dowson

# This file here is really only extra for experts (and because @lkapelevich asked for it)
function sample(x::Vector, w=ones(length(x))/length(x))
    @assert length(w) == length(x)
    @assert (sum(w) - 1.0) <= 1e-6
    r = rand(Float64)
    @inbounds for i in 1:length(x)
        r -= w[i]
        if r < 1e-6
            return x[i]
        end
    end
end
function sampleprob(x::Vector)
    @assert (sum(x) - 1.0) <= 1e-6
    r = rand(Float64)
    @inbounds for i in 1:length(x)
        r -= x[i]
        if r < 1e-6
            return i
        end
    end
end

abstract AbstractForwardSampler

samplemarkov(sampler::AbstractForwardSampler, m::SDDPModel, stage, lastmarkov) = error("""
    You must define a samplemarkov method for your forward sampler $(typeof(sampler)).
""")

samplescenario(sampler::AbstractForwardSampler, scenarios::Vector{Scenario}, m::SDDPModel, stage, markov) = error("""
    You must define a samplescenario method for your forward sampler $(typeof(sampler)).
""")

samplepricescenario(sampler::AbstractForwardSampler, pricescenarios::Vector{PriceScenario}, m::SDDPModel, stage, markov) = error("""
    You must define a samplepricescenario method for your forward sampler $(typeof(sampler)).
""")

immutable MonteCarloSampler <: AbstractForwardSampler
    # transition[stage][i, j]
    transition::Vector{Arrray{Float64, 2}}
    # initialmarkovprobability[i]
    initialmarkovprobability::Vector{Float64}

    # scenarioprobability[stage][markov][scenario]
    scenarioprobabilitiy::Vector{Vector{Vector{Float64}}}

    # pricescenarioprobability[stage][markov][scenario]
    pricescenarioprobabilitiy::Vector{Vector{Vector{Float64}}}
end

function samplemarkov(sampler::MonteCarloSampler, m::SDDPModel, stage, lastmarkov)
    if stage == 1
        sampleprob(sampler.initialmarkovprobability)
    else
        sampleprob(sampler.transition[stage][lastmarkov, :])
    end
end

function samplescenario(sampler::MonteCarloSampler, scenarios::Vector{Scenario}, m::SDDPModel, stage, markov)
    idx = sampleprob(sample.scenarioprobability[stage][markov])
    return scenarios[idx]
end

function samplescenario(sampler::MonteCarloSampler, scenarios::Vector{PriceScenario}, m::SDDPModel, stage, markov)
    idx = sampleprob(sample.pricescenarioprobability[stage][markov])
    return scenarios[idx]
end

"""
    Sample the forward pass with the a uniform probability across all outcomes.
"""
immutable UniformSampler <: AbstractForwardSampler end

samplemarkov(sampler::UniformSampler, m::SDDPModel, stage, lastmarkov) = sample(1:m.markovstates[stage])

samplescenario(sampler::UniformSampler, scenarios::Vector{Scenario}, m::SDDPModel, stage, markov) = sample(scenarios)

samplescenario(sampler::UniformSampler, scenarios::Vector{PriceScenario}, m::SDDPModel, stage, markov) = sample(scenarios)
