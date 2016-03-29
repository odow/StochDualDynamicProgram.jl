# TODO
#  - at the moment we assume uniform scenario probability in each markov state
#  - cut selection [de Matos, Philpott, Finardi (2015). Improving the performance of stochastic dual dynamic programming]
#  - cut addition. why add a cut twice?

module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp
using Formatting
using Distributions

export SDDPModel,
    @defStateVar, @defValueToGo, @addScenarioConstraint, @setStageProfit,
    simulate, load_cuts!

include("macros.jl")
include("SDDPModel.jl")
include("SDDPalgorithm.jl")

"""
Solve the model using the SDDP algorithm.

Inputs:
m                  - the SDDP model object
simulation_passes  - the number of realisations to conduct when testing for convergence
maximum_iterations - the maximum number of iterations (cutting passes, convergence testing) to complete before termination
log_frequency      - simulate bound (using n=simulation_passes) every [log_frequency] iterations and output to user
beta_quantile      - CVar quantile for nested risk aversion
risk_lambda        - Weighting on convex combination of Expectation and CVar
    risk_lambda * Expectation + (1 - risk_lambda) * CVar
"""
function JuMP.solve{M,N,S,T}(m::SDDPModel{M,N,S,T}; simulation_passes=1, log_frequency=1, maximum_iterations=1, beta_quantile=1, risk_lambda=1)
    print_stats_header()

    # Set risk aversion parameters
    m.beta_quantile = beta_quantile
    m.risk_lambda = risk_lambda

    for i =1:maximum_iterations
        # Cutting passes
        backward_pass!(m)

        if mod(i, log_frequency) == 0
            # Simulate
            forward_pass!(m, simulation_passes)
            print_stats(m)
        end
    end

end

function print_stats(m::SDDPModel)
    printfmt("({1:8.2f}, {2:8.2f}) | {3:8.2f} | {4:3.2f}\n", m.confidence_interval[1], m.confidence_interval[2], m.valid_bound, 100*rtol(m))
end
function print_stats_header()
    printfmt("{1:22s} | {2:10s} | {3:6s}\n", "Expected Objective", "Valid Bound", "% Gap")
end

end
