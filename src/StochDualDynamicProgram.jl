#  Copyright 2017, Oscar Dowson

module StochDualDynamicProgram

export SDDPModel,
    @state, @states, @scenario, @scenarioconstraints, @stageobjective,
    objectivescenario!,
    DiscreteDistribution,
    Expectation, NestedCVaR

include("type_definitions.jl")
include("macro_definitions.jl")
include("function_definitions.jl")
include("cut_selection.jl")
include("risk_measures.jl")

end
