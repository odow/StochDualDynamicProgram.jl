#  Copyright 2016, Oscar Dowson

# Set the risk measure for all the subproblems
function setriskmeasure!(m::SDDPModel, riskmeasure::AbstractRiskMeasure)
    for sp in m.stage_problems
        stagedata(sp).beta_quantile = riskmeasure.beta
        stagedata(sp).lambda_weight = riskmeasure.lambda
    end
end

"""
    reweightscenarios!(m, t, i, beta, lambda)

This function recomputes the scenario weightings for stage `t` markov state `i`.
    ⭔ρ[z] = λ E[z] + (1-λ) CVaR[z]

Inputs:
    m       the SDDDP model object
    t       the current stage
    i       the current markov state
    beta    the CVaR β quantile
    lambda  the λ weight on expectation
"""
function reweightscenarios!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, obj::Vector, t::Int, i::Int, beta::Float64, lambda::Float64)
    if abs(lambda - 1) > 1e-5     # We are risk averse
        for s=1:S
            for j=1:M
                stagedata(m, t, i).weightings_matrix[j,s] = transitionprobability(m, t, i, j)*m.scenario_probability[s]
            end
        end
        # Reweighted probabilities
        nestedcvar!(stagedata(m, t, i).weightings_matrix, obj, beta, lambda, X)
        @assert abs(sum(stagedata(m, t, i).weightings_matrix) - 1) < 1e-5 # Check sanity
    end
    return
end
function reweightscenarios!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t, i, beta, lambda)
    obj = zeros(M*S)
    idx = 1
    for s=1:S
        for j=1:M
            obj[idx] = stagedata(m, t+1, j).objective_values[s]
            idx += 1
        end
    end
    reweightscenarios!(m, obj, t, i, beta, lambda)
end

"""
    nestedcvar!(prob, x, beta, lambda, maximisation)

This function computes the risk averse weightings for CVaR given the beta quantile

Inputs:
    x      objective values for scenarios
    p      probability density for scenarios
    ß      CVar beta quantile
    ismax  (true/false) problem is maximisation
"""
function nestedcvar!{T}(p::Array{T}, x::Vector{T},  beta::Float64, lambda::Float64, ismax::Bool)
    @assert length(p) == length(x)                  # sanity
    q = 0.                                          # Quantile collected so far
    for i in sortperm(x, rev=!ismax)                # For each scenario in order
        if q >=  beta                               # We have collected the beta quantile
            riskaverse  = 0.                        # riskaverse contribution is zero
        else
            riskaverse  = min(p[i], beta-q) / beta  # risk averse contribution
        end
        p[i] = lambda * p[i] + (1 - lambda) * riskaverse  # take the biggest proportion of the scenario possible
        q += riskaverse * beta                      # Update total quantile collected
    end
end
nestedcvar!{T}(p::Array{T},  x::Vector{T}, beta::Float64, lambda::Float64, ::Type{Max}) = nestedcvar!(p,x,beta,lambda,true)
nestedcvar!{T}(p::Array{T},  x::Vector{T}, beta::Float64, lambda::Float64, ::Type{Min}) = nestedcvar!(p,x,beta,lambda,false)
