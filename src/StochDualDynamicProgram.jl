#  Copyright 2016, Oscar Dowson

module StochDualDynamicProgram

export SDDPModel,
    @state, @states, @scenario, @scenarioconstraints, @stageobjective,
    objectivescenario!,
    setriskmeasure!,
    DiscreteDistribution,
    Expectation, NestedCVaR

include("type_definitions.jl")
include("macro_definitions.jl")
include("function_definitions.jl")
include("cut_definitions.jl")

end
