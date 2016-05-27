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
function reweightscenarios!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t::Int, i::Int, beta::Float64, lambda::Float64)
    if abs(lambda - 1) > 1e-5     # We are risk averse
        P = zeros(S*M)            # Initialise storage
        idx=1
        for j=1:M
            for s=1:S
                P[idx] = transitionprobability(m, t, i, j)*m.scenario_probability[s]
                idx+=1
            end
        end
        # Reweighted probabilities
        nestedcvar!(P,
            vcat([stagedata(sp).objective_values for sp in m.stage_problems[t+1, :]]...),
            beta, lambda, X)
        idx=1
        for j=1:M
            for s=1:S
                stagedata(m, t, i).weightings_matrix[j, s] = P[idx]
                idx+=1
            end
        end
        @assert abs(sum(stagedata(m, t, i).weightings_matrix) - 1) < 1e-5 # Check sanity
    end
    return
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
function nestedcvar!{T}(p::Vector{T}, x::Vector{T},  beta::Float64, lambda::Float64, ismax::Bool)
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
nestedcvar!{T}(p::Vector{T},  x::Vector{T}, beta::Float64, lambda::Float64, ::Type{Max}) = nestedcvar!(p,x,beta,lambda,true)
nestedcvar!{T}(p::Vector{T},  x::Vector{T}, beta::Float64, lambda::Float64, ::Type{Min}) = nestedcvar!(p,x,beta,lambda,false)
