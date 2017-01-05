# Copyright Oscar Dowson, 2016

Cut(intercept, coefficients::Vector) = Cut(intercept, tuple(coefficients...))

function _dot{N}(x::Tuple{Vararg{N, Float64}}, y::Tuple{Vararg{N, Float64}})
    z = 0.0
    @inbounds for i=1:N
        z += x[i] * y[i]
    end
    z
end
evaluate(cut::Cut{N}, state::Tuple{Vararg{N, Float64}}) = cut.intercept + _dot(cut.coefficients, state)

function addcut!{N}(sense::Sense, cs::CutStorage{N}, csinner::CutStorageInnter{N}, cut::Cut{N})
    push!(csinner.cuts, cut)
    push!(csinner.states_dominant, 0)
    for i in 1:length(cs.statesvisited)
        val = evaluate(cut, cs.statesvisited[i])
        if dominates(val, csinner.best_bound[i])
            csinner.best_bound[i] = val
            csinner.best_cut_index[i] = length(csinner.cuts)
            csinner.states_dominant[end] += 1
        end
    end
end

# function rebuild!(::DeMatosCutSelection, m::SDDPModel)
# end

# TODO: need to be able to run cut selection on a single rib

"""
    This function assembles a new cut using the following inputs
    + oracle::CutOracle - used to dispatch
    + sense::Sense - either ::Maximum or ::Minimum
    + x::Vector{Vector{Float64}} - a vector of vector of state values for each scenario
    + pi::Vector{Vector{Float64}} - a vector of vector of dual values for each scenario
    + theta::Vector{Float64} - a vector of objective values for each scenario
    + prob::Vector{Float64} - the probability support of the scenarios. Should sum to one
    + stage::Int - the index of the stage
    + markov::Int - the index of the markov state
"""
cutgenerator(oracle::CutOracle, sense::AbstactSense, x, pi, theta, prob, stage, markov) =
    cutgenerator(oracle::CutOracle, sense::AbstactSense, x, pi, theta, prob)

cutgenerator(oracle::CutOracle, sense::AbstactSense, x, pi, theta, prob) = error("""
    You need to overload a `cutgenerator` method for the oracle of type $(typeof(oracle)).
    This could be the method including the stage and markov index
        cutgenerator(oracle::CutOracle, sense::AbstactSense, x, pi, theta, prob, stage, markov)
    or
        cutgenerator(oracle::CutOracle, sense::AbstactSense, x, pi, theta, prob)
""")

function cutgenerator(ex::Expectation, sense::AbstactSense, x, pi, theta, prob)
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
function cutgenerator(cvar::NestedCVaR, sense::AbstractSense, x, pi, theta, prob)
    @assert length(pi) == length(theta) == length(prob)
    if length(cvar.storage) < length(prob)
        append!(cvar.storage, zeros(length(prob) - length(cvar.storage)))
    end
    calculateCVaRprobabilities!(cvar.storage, sense, prob, x, cvar.lambda, cvar.beta)

    cutgenerator(expectation, x, pi, theta, cvar.storage[1:length(prob)])
end

_sortperm(::Maximum, x) = sortperm(x, rev=false)
_sortperm(::Minimum, x) = sortperm(x, rev=true)
function calculateCVaRprobabilities!(newprob, sense, oldprob, theta, lambda, beta::Float64)
    @assert length(newprob) >= length(oldprob)
    quantile_collected = 0.0
    cvarprob = 0.0
    cache = (1 - lambda) / beta             # cache this to save some operations
    @inbounds for i in _sortperm(sense, theta)        # For each scenario in order
        newprob[i] = lambda * oldprob[i]    # expectation contribution
        if quantile_collected <  beta       # We haven't collected the beta quantile
            cvarprob = min(oldprob[i], beta-quantile_collected)
            newprob[i] += cache * cvarprob  # risk averse contribution
            quantile_collected += cvarprob  # Update total quantile collected
        end
    end
end
