#  Copyright 2017, Oscar Dowson

module StochDualDynamicProgram

using JuMP

export SDDPModel,
    @state, @states, @scenario, @scenarioconstraints,
    stageobjective!, priceprocess!,
    # risk measures
    Expectation, NestedCVaR,
    # cut oracles
    DefaultCutOracle, DeMatosCutOracle
#     # forward samplers
#     DefaultSampler, UniformSampler

include("types.jl")
include("utilities.jl")

include("cutoracles.jl")
include("riskmeasures.jl")

include("scenario.jl")
include("statevariable.jl")
include("price.jl")

include("model.jl")

end
