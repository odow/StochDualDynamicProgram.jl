isdefined(Base, :__precompile__) && __precompile__()

module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp
using Formatting
using Distributions, StatsBase

import Base.dot

export SDDPModel,
    @state, @scenarioconstraint, @scenarioconstraints, @stageprofit,
    @visualise,
    simulate, loadcuts!,
    LevelOne, Deterministic, NoSelection,
    ConvergenceTest, BackwardPass, Parallel,
    NestedCVar,
    Convergence,
    NoRegularisation, LinearRegularisation, QuadraticRegularisation

include("types.jl")
include("macros.jl")
include("model.jl")
include("forwardpass.jl")
include("backwardpass.jl")
include("cut_selection.jl")
include("print.jl")
include("algorithm.jl")
include("simulate.jl")

end
