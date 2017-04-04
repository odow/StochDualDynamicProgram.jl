#  Copyright 2017, Oscar Dowson

using StochDualDynamicProgram, JuMP, Base.Test

const EXAMPLESDIR = joinpath(dirname(dirname(@__FILE__)), "examples")

@testset "SDDPModel" begin
    m = SDDPModel(
        stages=3, markov_states=2, transition=[0.2 0.8;0.8 0.2], value_to_go_bound=1000
        ) do sp, stage, markov
        nothing
    end
    @test size(m.stage_problems) == (3,2)
    @test StochDualDynamicProgram.transitionprobability(m, 1, 1, 1) == 0.2

    sp = StochDualDynamicProgram.StageProblem()
    @test StochDualDynamicProgram.issubproblem(sp) == true
    @test StochDualDynamicProgram.issubproblem(m) == false

    m = SDDPModel(markov_states=2, value_to_go_bound=1) do sp, stage, markov
        nothing
    end
    @test StochDualDynamicProgram.transitionprobability(m, 1, 1, 1) == 0.5
end

@testset "@state" begin
    m = StochDualDynamicProgram.StageProblem()
    @state(m, 0 <= x <= 3, x0=2.5)
    @test_approx_eq m.colLower [0, -Inf]
    @test_approx_eq m.colUpper [3, Inf]
    @test length(StochDualDynamicProgram.stagedata(m).dual_constraints) == 1
    @test getlowerbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) == 2.5
    @test getupperbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) == 2.5

    m = StochDualDynamicProgram.StageProblem()
    @state(m, y <= 1, y0=1.0)
    @test_approx_eq m.colLower [-Inf, -Inf]
    @test_approx_eq m.colUpper [1, Inf]
    @test length(StochDualDynamicProgram.stagedata(m).dual_constraints) == 1
    @test getlowerbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) == 1.0
    @test getupperbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) == 1.0

    m = StochDualDynamicProgram.StageProblem()
    @state(m, z, z0=1.0)
    @test_approx_eq m.colLower [-Inf, -Inf]
    @test_approx_eq m.colUpper [Inf, Inf]
    @test length(StochDualDynamicProgram.stagedata(m).dual_constraints) == 1
    @test getlowerbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) == 1.0
    @test getupperbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) == 1.0

    m = StochDualDynamicProgram.StageProblem()
    T = [:a, :b]
    rhs = Dict{Symbol, Float64}(:a=>1., :b=>2.)
    @state(m, x[t=T] >= 0., x0=rhs[t])
    @test_approx_eq m.colLower [0., 0., -Inf, -Inf]
    @test_approx_eq m.colUpper [Inf, Inf, Inf, Inf]
    @test length(StochDualDynamicProgram.stagedata(m).dual_constraints) == 2
    @test getlowerbound(StochDualDynamicProgram.stagedata(m).dual_constraints[1]) == rhs[:a]
    @test getupperbound(StochDualDynamicProgram.stagedata(m).dual_constraints[2]) == rhs[:b]

end

@testset "Hydro Example" begin
    include(joinpath(EXAMPLESDIR, "hydro.jl"))
    @test_approx_eq_eps mean(results[:Objective]) 904 20

    include(joinpath(EXAMPLESDIR, "hydro2.jl"))
    @test_approx_eq_eps mean(results[:Objective]) -1450 20
end

@testset "Newsvendor Example" begin
    include(joinpath(EXAMPLESDIR, "newsvendor.jl"))
    @test_approx_eq_eps mean(results[:Objective]) 97.5 1

    include(joinpath(EXAMPLESDIR, "newsvendor2.jl"))
    @test_approx_eq_eps mean(results[:Objective]) 99.36 1
end

@testset "Visualisation" begin
    include(joinpath(EXAMPLESDIR, "visualisation.jl"))
end

@testset "Parallelisation" begin
    include(joinpath(EXAMPLESDIR, "parallel_example.jl"))
    @test_approx_eq_eps mean(results[:Objective]) 906.35 5

    include(joinpath(EXAMPLESDIR, "serial_comparison.jl"))
    @test_approx_eq_eps mean(results[:Objective]) 906.35 5
end

@testset "Complete Example" begin
    include(joinpath(EXAMPLESDIR, "complete_example.jl"))
    @test_approx_eq_eps mean(results[:Objective]) 910 5
end

@testset "SOSII Example" begin
    include(joinpath(EXAMPLESDIR, "sos_hydro.jl"))
    @test_approx_eq_eps mean(results[:Objective]) -920 20
end
