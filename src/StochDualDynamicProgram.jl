#  Copyright 2017, Oscar Dowson

module StochDualDynamicProgram

using JuMP

export SDDPModel,
    @state, @states, @scenario, @scenarioconstraints, @stageobjective,
    objectivescenario!,
    DiscreteDistribution,
    Expectation, NestedCVaR

include("type_definitions.jl")
include("utilities.jl")
include("cut_selection.jl")
include("risk_measures.jl")
include("macro_definitions.jl")
include("function_definitions.jl")
include("model.jl")


end
