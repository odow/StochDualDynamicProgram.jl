#  Copyright 2017, Oscar Dowson

initialise(t::AbstractCutOracle, m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int) = error("""
    You must define an initialisation method for your cut oracle.
""")

"""
    storecut!(oracle::AbstactCutOracle, stage::Int, markovstate::Int, pricestate::Int, cut)

    This function adds the cut to the CutOracle.
"""
storecut!(oracle::AbstractCutOracle, m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut) = error("""
    You must define the function
        storecut!(oracle::$(typeof(oracle)), m::SDDPModel, stage::Int, markovstate::Int, pricestate::Int, cut)
    that is overloaded for your oracle of type $(typeof(oracle)).
""")

"""
    validcuts(oracle::AbstactCutOracle)

    This function returns an iterable list of all the valid cuts contained within the oracle.
"""
validcuts(oracle::AbstractCutOracle) = error("""
You must define the function validcuts(oracle::$(typeof(oracle))) that is overloaded for your
    oracle of type $(typeof(oracle)).
""")

include("cutoracles/default.jl")
include("cutoracles/dematos.jl")
