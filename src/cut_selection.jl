#  Copyright 2017, Oscar Dowson

"""
    addcut!(oracle::CutOracle, stage::Int, markovstate::Int, pricestate::Int, cut)

    This function adds the cut to the CutOracle.
"""
addcut!(oracle::CutOracle, stage::Int, markovstate::Int, pricestate::Int, cut) = error("""
    You must define the function
        addcut!(oracle::$(typeof(oracle)), stage::Int, markovstate::Int, pricestate::Int, cut)
    that is overloaded for your oracle of type $(typeof(oracle)).
""")

"""
    An optional method to overload if the oracle needs it
"""
addvisitedstate!(oracle::CutOracle, stage::Int, markovstate::Int, pricestate::Int, state::StateTuple{N}) = nothing

"""
    validcuts(oracle::CutOracle)

    This function returns an iterable list of all the valid cuts contained within the oracle.
"""
validcuts(oracle::CutOracle, stage::Int, markovstate::Int, pricestate::Int) = error("""
You must define the function validcuts(oracle::$(typeof(oracle))) that is overloaded for your
    oracle of type $(typeof(oracle)).
""")

"""
    A default cut storage oracle that simply remembers all the cuts
"""
immutable DefaultCutOracle{N} <: CutOracle
    cuts::Nested3Vector{Vector{Cut{N}}}
end
# DefaultCutOracle()

addcut!{N|}(oracle::DefaultCutOracle{N}, stage::Int, markovstate::Int, pricestate::Int, cut::Cut{N}) = push!(oracle.cuts[stage][markovstate][pricestate], cut)
validcuts{N}(oracle::DefaultCutOracle{N}, stage::Int, markovstate::Int, pricestate::Int) = oracle.cuts[stage][markovstate][pricestate]

immutable DeMatosCutOracle{N, S} <: CutOracle
    # statesvisited[state] = vector of state tuples
    statesvisited::Vector{Vector{StateTuple{N}}}
    # storage[stage][markov state][rib] = cut storage
    storage::Nested3Vector{DeMatosCutStorage{N, S}}
end
"""
    Vitor de Matos Level One Cut Selection
"""
immutable DeMatosCutStorage{N, S}
    cuts::Vector{Cut{N}}
    states_dominant::Vector{Int} # states_dominant[i] = number of states in statesvisited that cut i is dominant
    best_cut_index::Vector{Int} # best_cut_index[i] = index of cut in cuts that is the dominant cut at states_dominant[i]
    best_bound::Vector{Float64} # best_bound[i] = best objective bound at statesvisited[i]
end

function _dot{N, T}(x::Tuple{Vararg{N, T}}, y::Tuple{Vararg{N, T}})
    z = zero(T)
    @inbounds for i=1:N
        z += x[i] * y[i]
    end
    z
end
evaluate(cut::Cut{N}, state::StateTuple{N}) = cut.intercept + _dot(cut.coefficients, state)
dominates(::Maximisation, x, y) = x < y
dominates(::Minimisation, x, y) = x > y
bestval(::Maximisation) = Inf
bestval(::Minimisation) = -Inf

function _addcut!{N, S}(dematos::DeMatosCutStorage{N, S}, statesvisited::Vector{StateTuple{N}}, cut::Cut{N})
    push!(dematos.cuts, cut)
    push!(dematos.states_dominant, 0)
    for i in 1:length(statesvisited)
        val = evaluate(cut, statesvisited[i])
        if dominates(S, val, dematos.best_bound[i])
            dematos.best_bound[i] = val
            dematos.best_cut_index[i] = length(dematos.cuts)
            dematos.states_dominant[end] += 1
        end
    end
end
addcut!{N, S}(dematos::DeMatosCutOracle{N, S}, stage::Int, markovstate::Int, pricestate::Int, cut::Cut{N}) = _addcut!(dematos.storage[stage][markovstate][pricestate], dematos.statesvisited[stage], cut)

function _addvisitedstate!{N, S}(dematos::DeMatosCutStorage{N, S}, statesvisited::Vector{Vector{Vector{DeMatosCutStorage{N, S}}}, state::StateTuple{N})
    push!(statesvisited, state)
    dematos.best_bound[end] = bestval(S)
    dematos.best_cut_index[end] = 0
    for i in 1:length(dematos.cuts)
        val = evaluate(dematos.cuts[i], statesvisited[end])
        if dominates(S, val, dematos.best_bound[end])
            dematos.best_bound[end]     = val
            dematos.best_cut_index[end] = i
        end
    end
    dematos.states_dominant[dematos.best_cut_index[end]] += 1
end
addvisitedstate!{N, S}(dematos::DeMatosCutOracle{N, S}, stage::Int, markovstate::Int, pricestate::Int, state::StateTuple{N}) = addvisitedstate!(dematos.storage[stage][markovstate][pricestate], dematos.statesvisited[stage], state)

Base.start{N, S}(dematos::DeMatosCutStorage{N, S}) = 1
Base.done{N, S}(dematos::DeMatosCutStorage{N, S}, state) = (state == -1)
function Base.next{N, S}(dematos::DeMatosCutStorage{N, S}, state)
    @inline for i=state:length(dematos.cuts)
        if dematos.states_dominant[i] > 0
            return (dematos.cuts[i], i+1)
        end
    end
    return (Cut(N), -1)
end
validcuts{N, S}(dematos::DeMatosCutSelection{N, S}, stage::Int, markovstate::Int, pricestate::Int) = dematos.storage[stage][markovstate][pricestate]
