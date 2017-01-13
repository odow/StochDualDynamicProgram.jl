immutable NonRiskMeasure end

@testset "Risk Measures" begin
    m = SDDPModel(SDDP.Maximisation, ()->())
    @test_throws Exception SDDP.cutgenerator(NonRiskMeasure(), m, zeros(2), Float64[], Float64[], Float64[])

    x     = [1, 1]
    prob  = [0.25, 0.25, 0.25, 0.25]
    theta = [3.0, 2.0, 1.0, 4.0]
    pi    = [ [0.5, 0.0], [0.0, 0.5], [0.5, 0.5], [1.0, 2.0] ]

    @testset "Expectation" begin
        cut = SDDP.cutgenerator(Expectation(), m, x, pi, theta, prob, 1, 1)
        @test cut.intercept == 1.25
        @test cut.coefficients == [0.5, 0.75]
    end

    @testset "NestedCV@R" begin
        cut1 = SDDP.cutgenerator(NestedCVaR(lambda=1.0, beta=0.5), m, x, pi, theta, prob, 1, 1)
        @test cut1.intercept == 1.25
        @test cut1.coefficients == [0.5, 0.75]

        cut2 = SDDP.cutgenerator(NestedCVaR(lambda=0, beta=0.5), m, x, pi, theta, prob, 1, 1)
        @test cut2.intercept == 0.75
        @test cut2.coefficients == [0.25, 0.5]

        cut3 = SDDP.cutgenerator(NestedCVaR(lambda=0.5, beta=0.5), m, x, pi, theta, prob, 1, 1)
        @test cut3.intercept == 1.0
        @test cut3.coefficients == [0.375, 0.625]
    end

end
