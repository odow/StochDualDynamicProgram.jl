using StochDualDynamicProgram, JuMP, FactCheck

const EXAMPLESDIR = joinpath(dirname(dirname(@__FILE__)), "examples")

facts("SDDPModel") do
    m = SDDPModel(
        stages=3, markov_states=2, transition=[0.2 0.8;0.8 0.2], value_to_go_bound=1000
        ) do sp, stage, markov
        nothing
    end
    @fact size(m.stage_problems) --> (3,2)
    @fact StochDualDynamicProgram.transitionprobability(m, 1, 1, 1) --> 0.2

    sp = StochDualDynamicProgram.StageProblem()
    @fact StochDualDynamicProgram.issubproblem(sp) --> true
    @fact StochDualDynamicProgram.issubproblem(m) --> false

    @fact_throws SDDPModel(sense=:bogus)

    m = SDDPModel(markov_states=2, value_to_go_bound=1) do sp, stage, markov
        nothing
    end
    @fact StochDualDynamicProgram.transitionprobability(m, 1, 1, 1) --> 0.5
end

facts("@state") do
    m = StochDualDynamicProgram.StageProblem()
    @state(m, 0 <= x <= 3, x0=2.5)
    @fact m.colLower --> roughly([0, -Inf], 1e-4)
    @fact m.colUpper --> roughly([3, Inf], 1e-4)
    @fact length(StochDualDynamicProgram.stagedata(m).dual_constraints) --> 1
    @fact getlowerbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) --> 2.5
    @fact getupperbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) --> 2.5

    m = StochDualDynamicProgram.StageProblem()
    @state(m, y <= 1, y0=1.0)
    @fact m.colLower --> roughly([-Inf, -Inf], 1e-4)
    @fact m.colUpper --> roughly([1, Inf], 1e-4)
    @fact length(StochDualDynamicProgram.stagedata(m).dual_constraints) --> 1
    @fact getlowerbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) --> 1.0
    @fact getupperbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) --> 1.0

    m = StochDualDynamicProgram.StageProblem()
    @state(m, z, z0=1.0)
    @fact m.colLower --> roughly([-Inf, -Inf], 1e-4)
    @fact m.colUpper --> roughly([Inf, Inf], 1e-4)
    @fact length(StochDualDynamicProgram.stagedata(m).dual_constraints) --> 1
    @fact getlowerbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) --> 1.0
    @fact getupperbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) --> 1.0

    m = StochDualDynamicProgram.StageProblem()
    T = [:a, :b]
    rhs = Dict{Symbol, Float64}(:a=>1., :b=>2.)
    @state(m, x[t=T] >= 0., x0=rhs[t])
    @fact m.colLower --> roughly([0., 0., -Inf, -Inf], 1e-4)
    @fact m.colUpper --> roughly([Inf, Inf, Inf, Inf], 1e-4)
    @fact length(StochDualDynamicProgram.stagedata(m).dual_constraints) --> 2
    @fact getlowerbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) --> rhs[:a]
    @fact getupperbound(StochDualDynamicProgram.stagedata(m).dual_constraints[2]) --> rhs[:b]

end

facts("Hydro Example") do
    include(joinpath(EXAMPLESDIR, "hydro.jl"))
    @fact mean(results[:Objective]) --> roughly(904, 20)

    include(joinpath(EXAMPLESDIR, "hydro2.jl"))
    @fact mean(results[:Objective])--> roughly(-1450, 20)
end

facts("Newsvendor Example") do
    include(joinpath(EXAMPLESDIR, "newsvendor.jl"))
    @fact mean(results[:Objective]) --> roughly(97.5, 1)

    include(joinpath(EXAMPLESDIR, "newsvendor2.jl"))
    @fact mean(results[:Objective]) --> roughly(99.36, 1)
end

facts("Visualisation") do
    include(joinpath(EXAMPLESDIR, "visualisation.jl"))
end

facts("Parallelisation") do
    include(joinpath(EXAMPLESDIR, "parallel_example.jl"))
    @fact mean(results[:Objective]) --> roughly(906.35, 5)

    include(joinpath(EXAMPLESDIR, "serial_comparison.jl"))
    @fact mean(results[:Objective]) --> roughly(906.35, 5)
end

facts("Complete Example") do
    include(joinpath(EXAMPLESDIR, "complete_example.jl"))
    @fact mean(results[:Objective]) --> roughly(910, 5)
end

facts("SOSII Example") do
    include(joinpath(EXAMPLESDIR, "sos_hydro.jl"))
    @fact mean(results[:Objective]) --> roughly(-1450, 20)
end

FactCheck.exitstatus()
