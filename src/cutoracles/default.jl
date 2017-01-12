#  Copyright 2017, Oscar Dowson

"""
    A default cut storage oracle that simply remembers all the cuts
"""
immutable DefaultCutOracle <: AbstractCutOracle
    cuts::Vector{Cut}
end

DefaultCutOracle() = DefaultCutOracle(Cut[])

storecut!(oracle::DefaultCutOracle, m, stage, markovstate, pricestate, cut::Cut) = push!(oracle.cuts, cut)

validcuts(oracle::DefaultCutOracle) = oracle.cuts

initialise(::DefaultCutOracle, m, stage::Int, markovstate::Int, pricestate::Int) = DefaultCutOracle(Cut[])
