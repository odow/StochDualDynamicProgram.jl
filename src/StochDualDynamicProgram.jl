#  Copyright 2017, Oscar Dowson

module StochDualDynamicProgram

using JuMP

export SDDPModel,
    @state, @states, @scenario, @scenarioconstraints, @stageobjective,
    objectivescenario!,
    DiscreteDistribution,
    # risk measures
    Expectation, NestedCVaR,
    # cut oracles
    DefaultCutOracle, DeMatosCutOracle,
    # forward samplers
    DefaultSampler, UniformSampler

include("type_definitions.jl")
include("utilities.jl")
include("cut_selection.jl")
include("risk_measures.jl")
include("forward_sampler.jl")
include("macro_definitions.jl")
include("function_definitions.jl")
include("model.jl")


end
