using StochDualDynamicProgram, JuMP, FactCheck

facts("SDDPModel") do
    m = SDDPModel(
        stages=3,
        markov_states=2,
        transition=[0.2 0.8;0.8 0.2]
    )
    @fact size(m.stage_problems) --> (3,2)
    @fact StochDualDynamicProgram.get_transition(m, 1, 1, 1) --> 0.2

    sp = StochDualDynamicProgram.StageProblem()
    @fact StochDualDynamicProgram.is_sp(sp) --> true
    @fact StochDualDynamicProgram.is_sp(m) --> false
end

facts("@defStateVar") do
    m = StochDualDynamicProgram.StageProblem()
    @defStateVar(m, 0 <= x <= 3, x0==2.5)
    @fact :x in m.ext[:state_vars] --> true
    @fact m.colLower --> roughly([0, -Inf], 1e-4)
    @fact m.colUpper --> roughly([3, Inf], 1e-4)
    @fact haskey(m.ext[:duals], :x) --> true
    @fact getLower(m.ext[:duals][:x]) --> 2.5
    @fact getUpper(m.ext[:duals][:x]) --> 2.5

    m = StochDualDynamicProgram.StageProblem()
    @defStateVar(m, y <= 1, y0==1.0)
    @fact :y in m.ext[:state_vars] --> true
    @fact m.colLower --> roughly([-Inf, -Inf], 1e-4)
    @fact m.colUpper --> roughly([1, Inf], 1e-4)
    @fact haskey(m.ext[:duals], :y) --> true
    @fact getLower(m.ext[:duals][:y]) --> 1.0
    @fact getUpper(m.ext[:duals][:y]) --> 1.0

    m = StochDualDynamicProgram.StageProblem()
    @defStateVar(m, z, z0==1.0)
    @fact :z in m.ext[:state_vars] --> true
    @fact m.colLower --> roughly([-Inf, -Inf], 1e-4)
    @fact m.colUpper --> roughly([Inf, Inf], 1e-4)
    @fact haskey(m.ext[:duals], :z) --> true
    @fact getLower(m.ext[:duals][:z]) --> 1.0
    @fact getUpper(m.ext[:duals][:z]) --> 1.0
end

facts("@defValueToGo") do
    m = StochDualDynamicProgram.StageProblem()
    @defValueToGo(m, theta)
    @fact haskey(m.ext, :theta) --> true
    @fact typeof(m.ext[:theta]) --> Variable
    @fact m.colLower --> roughly([-Inf], 1e-4)
    @fact m.colUpper --> roughly([Inf], 1e-4)

    m = StochDualDynamicProgram.StageProblem()
    @defValueToGo(m, theta <= 10)
    @fact haskey(m.ext, :theta) --> true
    @fact typeof(m.ext[:theta]) --> Variable
    @fact m.colLower --> roughly([-Inf], 1e-4)
    @fact m.colUpper --> roughly([10], 1e-4)

    m = StochDualDynamicProgram.StageProblem()
    @defValueToGo(m, 0 <= theta <= 10)
    @fact haskey(m.ext, :theta) --> true
    @fact typeof(m.ext[:theta]) --> Variable
    @fact m.colLower --> roughly([0], 1e-4)
    @fact m.colUpper --> roughly([10], 1e-4)
end

facts("Hydro Example") do
    context("Version One") do
        include("../examples/hydro.jl")
        results = solve_hydro()
        @fact Float64[mean(results[:upper_reservoir][t]) for t=1:3] --> roughly([122.99, 66.83, 0.0], 0.1)
        @fact Float64[mean(results[:lower_reservoir][t]) for t=1:3] --> roughly([200.0, 75.02, 0.0], 0.1)
    end

    context("Version Two") do
        include("../examples/hydro2.jl")
        results = solve_hydro2()
        @fact Float64[mean(results[:upper_reservoir][t]) for t=1:3] --> roughly([122.99, 66.83, 0.0], 0.1)
        @fact Float64[mean(results[:lower_reservoir][t]) for t=1:3] --> roughly([200.0, 75.02, 0.0], 0.1)
    end
end
