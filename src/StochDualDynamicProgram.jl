# TODO at the moment we assume uniform scenario probability in each markov state

module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp
using Formatting
using Distributions

export SDDPModel,
    @defStateVar, @defValueToGo, @addScenarioConstraint, @setStageProfit,
    simulate

include("macros.jl")
include("SDDPModel.jl")
include("SDDPalgorithm.jl")

"""
Solve the model using the SDDP algorithm.

Inputs:
m               - the SDDP model object
forward_passes  - the number of realisations to conduct when testing for convergence
backward_passes - the number of cuttting passes to conduct before testing for convergence
max_iterations  - the maximum number of iterations (cutting passes, convergence testing) to complete before termination
beta_quantile   - CVar quantile for nested risk aversion
risk_lambda     - Weighting on convex combination of Expectation and CVar
    risk_lambda * Expectation + (1 - risk_lambda) * CVar
"""
function JuMP.solve{M,N,S,T}(m::SDDPModel{M,N,S,T}; forward_passes=1, backward_passes=1, max_iterations=1000, beta_quantile=1, risk_lambda=1)
    print_stats_header()

    # Set risk aversion parameters
    m.beta_quantile = beta_quantile
    m.risk_lambda = risk_lambda

    i=0
    while i < max_iterations
        # Cutting passes
        backward_pass!(m, backward_passes)

        # Simulate
        forward_pass!(m, forward_passes)

        print_stats(m)

        i += 1
    end

end

function print_stats(m::SDDPModel)
    printfmt("({1:8.2f}, {2:8.2f}) | {3:8.2f} | {4:3.2f}\n", m.confidence_interval[1], m.confidence_interval[2], m.valid_bound, 100*rtol(m))
end
function print_stats_header()
    printfmt("{1:22s} | {2:10s} | {3:6s}\n", "Expected Objective", "Valid Bound", "% Gap")
end

end
