# TODO
#  - at the moment we assume uniform scenario probability in each markov state
#  - cut selection [de Matos, Philpott, Finardi (2015). Improving the performance of stochastic dual dynamic programming]

module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp
using Formatting
using Distributions, StatsBase

export SDDPModel,
    @defStateVar, @defValueToGo, @addScenarioConstraint, @setStageProfit,
    simulate, load_cuts!

include("macros.jl")
include("de_matos_types.jl")
include("SDDPModel.jl")
include("de_matos_functions.jl")
include("SDDPalgorithm.jl")
include("parallel_backwards.jl")
include("parallel_forewards.jl")


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
check_duplicate_cuts - Checks cut before adding to see if duplicate
"""
# function JuMP.solve(m::SDDPModel; simulation_passes=1, log_frequency=1, maximum_iterations=1, beta_quantile=1, risk_lambda=1,  cut_selection_frequency=0)
#     print_stats_header()
#
#     cut_selection_frequency > 0 && initialise_cut_selection!(m)
#     # Set risk aversion parameters
#     m.beta_quantile = beta_quantile
#     m.risk_lambda = risk_lambda
#     # try
#         for i =1:maximum_iterations
#             # Cutting passes
#             backward_pass!(m, cut_selection_frequency > 0)
#
#             if mod(i, log_frequency) == 0
#                 # Simulate
#                 (_flag, n) = forward_pass!(m, simulation_passes)
#                 print_stats(m, n)
#                 !_flag && return
#             end
#             if cut_selection_frequency > 0 && mod(i, cut_selection_frequency) == 0
#                 rebuild_stageproblems!(m)
#             end
#         end
#     # catch InterruptException
#         # warn("Terminating early")
#     # end
#     return
# end
function JuMP.solve(m::SDDPModel; simulation_passes=1, log_frequency=1, maximum_iterations=1, beta_quantile=1, risk_lambda=1,  cut_selection_frequency=0)
    length(workers()) > 1 && initialise_workers!(m)

    print_stats_header()

    initialise_cut_selection!(m)
    # Set risk aversion parameters
    m.beta_quantile = beta_quantile
    m.risk_lambda = risk_lambda
    # try
    iterations = 0
    while iterations < maximum_iterations
        # Cutting passes
        n = min(max(1,log_frequency), maximum_iterations-iterations)
        parallel_backward_pass!(m, n)

        iterations += n

        # if mod(i, log_frequency) == 0
            # Simulate
        (_flag, n) = parallel_forward_pass!(m, simulation_passes)
        print_stats(m, n)
        !_flag && return
        # end
    end
    # catch InterruptException
        # warn("Terminating early")
    # end
    return
end

function print_stats(m::SDDPModel, n)
    printfmt("({1:8.2f}, {2:8.2f}) | {3:8.2f} | {4:3.2f} | {5:10d}\n", m.confidence_interval[1], m.confidence_interval[2], m.valid_bound, 100*rtol(m), round(Int, n))
end
function print_stats_header()
    printfmt("{1:22s} | {2:10s} | {3:6s} | {4:10s}\n", "Expected Objective", "Valid Bound", "% Gap", "# Samples")
end

end
