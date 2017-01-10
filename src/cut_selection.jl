#  Copyright 2017, Oscar Dowson
initialise(t::CutOracle, sense::Sense, stage::Int, markovstate::Int, pricestate::Int) = error("""
    You must define an initialisation method for your cut oracle.
""")

"""
    storecut!(oracle::CutOracle, stage::Int, markovstate::Int, pricestate::Int, cut)

    This function adds the cut to the CutOracle.
"""
storecut!{M<:SDDPModel}(oracle::CutOracle, m::M, stage::Int, markovstate::Int, pricestate::Int, cut) = error("""
    You must define the function
        storecut!(oracle::$(typeof(oracle)), m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut)
    that is overloaded for your oracle of type $(typeof(oracle)).
""")

"""
    validcuts(oracle::CutOracle)

    This function returns an iterable list of all the valid cuts contained within the oracle.
"""
validcuts(oracle::CutOracle) = error("""
You must define the function validcuts(oracle::$(typeof(oracle))) that is overloaded for your
    oracle of type $(typeof(oracle)).
""")

"""
    A default cut storage oracle that simply remembers all the cuts
"""
immutable DefaultCutOracle <: CutOracle
    cuts::Vector{Cut}
end
DefaultCutOracle() = DefaultCutOracle(Cut[])
storecut!(oracle::DefaultCutOracle, m, stage, markovstate, pricestate, cut::Cut) = push!(oracle.cuts, cut)
validcuts(oracle::DefaultCutOracle) = oracle.cuts
initialise(::DefaultCutOracle, sense::Sense, stage::Int, markovstate::Int, pricestate::Int) = DefaultCutOracle(Cut[])


"""
    Vitor de Matos Level One Cut Selection
"""
immutable DeMatosCutOracle
    sense::Sense
    cuts::Vector{Cut}
    states_dominant::Vector{Int} # states_dominant[i] = number of states in statesvisited that cut i is dominant
    best_cut_index::Vector{Int} # best_cut_index[i] = index of cut in cuts that is the dominant cut at states_dominant[i]
    best_bound::Vector{Float64} # best_bound[i] = best objective bound at statesvisited[i]
end
DeMatosCutOracle() = DeMatosCutOracle(Minimisation, Cut[], Int[], Int[], Float64[])
initialise(::DeMatosCutOracle, sense, stage, markovstate, pricestate) = DeMatosCutOracle(sense, Cut[], Int[], Int[], Float64[])

evaluate(cut::Cut, state::Vector{Float64}) = cut.intercept + dot(cut.coefficients, state)
dominates(::Type{Maximisation}, x, y) = x < y
dominates(::Type{Minimisation}, x, y) = x > y
bestval(::Type{Maximisation}) = Inf
bestval(::Type{Minimisation}) = -Inf

function _storecut!(dematos::DeMatosCutOracle, statesvisited, cut)
    push!(dematos.cuts, cut)
    push!(dematos.states_dominant, 0)
    @inbounds for i in 1:length(dematos.best_bound)
        val = evaluate(cut, statesvisited[i])
        if dominates(dematos.sense, val, dematos.best_bound[i])
            dematos.best_bound[i] = val
            dematos.best_cut_index[i] = length(dematos.cuts)
            dematos.states_dominant[end] += 1
        end
    end
end
function _update_level_one_data!(dematos::DeMatosCutOracle, statesvisited)
    @inbounds for i in (length(dematos.best_bound) + 1):length(statesvisited)
        push!(dematos.best_bound, bestval(dematos.sense))
        push!(dematos.best_cut_index, 0)
        @inbounds for j in 1:length(dematos.cuts)
            val = evaluate(dematos.cuts[j], statesvisited[i])
            if dominates(dematos.sense, val, dematos.best_bound[j])
                dematos.best_bound[j]     = val
                dematos.best_cut_index[j] = i
            end
        end
        dematos.states_dominant[dematos.best_cut_index[i]] += 1
    end
end

function storecut!(dematos::DeMatosCutOracle, m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut::Cut)
    _storecut!(dematos, m.statesvisited[stage], cut)
    _update_level_one_data!(dematos, m.statesvisited[stage])
end

Base.start(dematos::DeMatosCutOracle) = 1
Base.done(dematos::DeMatosCutOracle, state) = (state == -1)
function Base.next(dematos::DeMatosCutOracle, state)
    @inbounds for i=state:length(dematos.cuts)
        if dematos.states_dominant[i] > 0
            return (dematos.cuts[i], i+1)
        end
    end
    return (Cut(), -1)
end
validcuts(dematos::DeMatosCutOracle) = dematos
