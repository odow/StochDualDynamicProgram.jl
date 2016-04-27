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
    simulate, load_cuts!,
    LevelOne, Deterministic

const CONVERGENCE_TERMINATION = :Convergenced
const ITERATION_TERMINATION = :MaximumIterations
const UNKNOWN_TERMINATION = :Unknown

abstract CutSelectionMethod
immutable LevelOne <: CutSelectionMethod end
immutable Deterministic <: CutSelectionMethod end
immutable NoSelection <: CutSelectionMethod end

include("macros.jl")
include("de_matos_types.jl")
include("SDDPModel.jl")
include("de_matos_functions.jl")
include("SDDPalgorithm.jl")
include("parallel.jl")
include("async.jl")


"""
Solve the model using the SDDP algorithm.

Inputs:
m::SDDPModel
    - the SDDP model object
simulation_passes::Int
    - the number of realisations to conduct when testing for convergence
maximum_iterations::Int
    - the maximum number of iterations (cutting passes, convergence testing) to complete before termination
convergence_test_frequency::Int
    - simulate bound (using n=simulation_passes) every [convergence_test_frequency] iterations and output to user
    - if [convergence_test_frequency] = 0, never test convergence. Terminate at [maximum_iterations]
beta_quantile::Float64 ∈ [0, 1]
    - CVar quantile for nested risk aversion
risk_lambda::Float64 ∈ [0, 1]
    - Weighting on convex combination of Expectation and CVar
        risk_lambda * Expectation + (1 - risk_lambda) * CVar
cut_selection_frequency::Int
    - Number of cutting passes to conduct before removing those cuts that are level one dominated
cuts_per_processor::Int
    - if [cuts_per_processor] > 0, [cuts_per_processor] cuts are computed on each available processors and combined to remove duplicates.
        in addition, convergence test (forward simulation) passes are divided up to each available processor and simulated in parallel.
    - if [cuts_per_processor] = 0, method runs in serial mode.
convergence_termination::Bool
    - If a convergence test is conducted with the bounds are found to have converged, and [convergence_termination] = true, method will terminate.
        If this is false, the method will terminate at [maximum_iterations]
"""
function JuMP.solve(m::SDDPModel; simulation_passes=1, convergence_test_frequency=1, maximum_iterations=1, beta_quantile=1, risk_lambda=1,  cut_selection_frequency=0, convergence_termination=false, cuts_per_processor=0, cut_selection_method=LevelOne())
    @assert beta_quantile >= 0 && beta_quantile <= 1
    @assert risk_lambda >= 0 && risk_lambda <= 1
    @assert convergence_test_frequency >= 0
    @assert maximum_iterations >= 0
    @assert cut_selection_frequency >= 0
    @assert cuts_per_processor >= 0

    parallel = (cuts_per_processor > 0)
    if parallel && length(workers()) < 2
        warn("Paralleisation requested (cuts_per_processor=$(cuts_per_processor)) but Julia is only running with a single processor. Start julia with `julia -p N` or use the `addprocs(N)` function to load N workers. Running in serial mode.")
        parallel = false
    end

    # Status
    status = UNKNOWN_TERMINATION

    # Initial output for user
    print_stats_header()

    # try
    if parallel
        status = parallel_solve(m, simulation_passes, convergence_test_frequency, maximum_iterations, beta_quantile, risk_lambda, cut_selection_frequency, convergence_termination, cuts_per_processor, cut_selection_method)
        # status = async_solve(m, simulation_passes, convergence_test_frequency, maximum_iterations, beta_quantile, risk_lambda, cut_selection_frequency, convergence_termination, cuts_per_processor)
    else
        status = serial_solve(m, simulation_passes, convergence_test_frequency, maximum_iterations, beta_quantile, risk_lambda, cut_selection_frequency, convergence_termination, cut_selection_method)
    end
    # catch InterruptException
        # warn("Terminating early")
    # end
    status
end

terminate(converged::Bool, do_test::Bool) = converged & do_test
terminate(converged::Bool, do_test::Bool, simulation_passes::Range) = converged
terminate(converged::Bool, do_test::Bool, simulation_passes) = terminate(converged, do_test)

function serial_solve{N1<:Real, N2<:Real}(m::SDDPModel, simulation_passes, convergence_test_frequency::Int, maximum_iterations::Int, beta_quantile::N1, risk_lambda::N2,  cut_selection_frequency::Int, convergence_test::Bool, cut_selection_method::CutSelectionMethod)
    # Initialise cut selection if appropriate
    if cut_selection_frequency > 0
        initialise_cut_selection!(m)
    end

    # Set risk aversion parameters
    m.beta_quantile, m.risk_lambda = beta_quantile, risk_lambda

    time_forwards = 0.
    time_backwards = 0.
    time_cut_selection = 0.
    simulations = 0
    for i=1:maximum_iterations
        # Cutting passes
        tic()
        backward_pass!(m, cut_selection_frequency > 0)
        time_backwards += toq()

        # Rebuild models if using Cut Selection
        if cut_selection_frequency > 0 && mod(i, cut_selection_frequency) == 0
            # @time deterministic_prune!(m)
            tic()
            rebuild_stageproblems!(m, cut_selection_method)
            time_cut_selection += toq()
        end

        # Test convergence if appropriate
        if convergence_test_frequency > 0 && mod(i, convergence_test_frequency) == 0
            tic()
            (is_converged, n) = forward_pass!(m, simulation_passes)
            time_forwards += toq()
            simulations += Int(n)

            # Output to user
            print_stats(m, i, time_backwards, simulations, time_forwards, time_cut_selection)

            # Terminate if converged and appropriate
            if terminate(is_converged, convergence_test, simulation_passes)
                return CONVERGENCE_TERMINATION
            end
        end
    end
    ITERATION_TERMINATION
end

function parallel_solve{N1<:Real, N2<:Real}(m::SDDPModel, simulation_passes::Int, convergence_test_frequency::Int, maximum_iterations::Int, beta_quantile::N1, risk_lambda::N2,  cut_selection_frequency::Int, convergence_test::Bool, cuts_per_processor::Int, cut_selection_method::CutSelectionMethod)
    # Intialise model on worker processors
    nworkers = initialise_workers!(m)

    # Initialise cut selection storage since we need it to communicate between processors
    initialise_cut_selection!(m)

    # Set risk aversion parameters
    m.beta_quantile, m.risk_lambda = beta_quantile, risk_lambda

    total_iterations = 0
    iterations_since_last_convergence_test = 0
    iterations_since_last_cut_selection = 0
    time_backwards = 0.
    time_forwards = 0.
    total_simulations = 0
    cut_selection_time = 0.
    while total_iterations < maximum_iterations
        # Update iteration counters
        total_iterations += cuts_per_processor * nworkers
        iterations_since_last_convergence_test += cuts_per_processor * nworkers
        iterations_since_last_cut_selection += cuts_per_processor * nworkers

        # Cutting passes
        tic()
        if cut_selection_frequency > 0 && iterations_since_last_cut_selection >= cut_selection_frequency
            cut_selection_time += parallel_backward_pass!(m, cuts_per_processor * nworkers, cut_selection_frequency, cut_selection_method)
            iterations_since_last_cut_selection = 0
        else
            parallel_backward_pass!(m, cuts_per_processor * nworkers, cut_selection_frequency, NoSelection())
        end
        time_backwards += toq()

        # Test convergence if appropriate
        if convergence_test_frequency > 0 && iterations_since_last_convergence_test >= convergence_test_frequency
            iterations_since_last_convergence_test = 0
            tic()
            (is_converged, n) = parallel_forward_pass!(m, simulation_passes)
            total_simulations += Int(n)
            time_forwards += toq()

            # Output to user
            print_stats(m, total_iterations, time_backwards, total_simulations, time_forwards, cut_selection_time)

            # Terminate if converged and appropriate
            if terminate(is_converged, convergence_test, simulation_passes)
                return CONVERGENCE_TERMINATION
            end
        end
    end
    ITERATION_TERMINATION
end

function print_stats(m::SDDPModel, iterations, back_time, simulations, forward_time, cut_time)
    printfmt("{1:>9s} {2:>9s} | {3:>9s} {4:>6.2f} | {5:6s} {6:6s} | {7:6s} {8:6s} | {9:5s}\n",
        humanize(m.confidence_interval[1], "7.3f"), humanize(m.confidence_interval[2], "7.3f"), humanize(m.valid_bound, "7.3f"), 100*rtol(m), humanize(iterations), humanize(back_time), humanize(simulations), humanize(forward_time), humanize(cut_time, "6.2f"))
end
function print_stats_header()
    printfmt("{1:38s} | {2:13s} | {3:13s} |   Cut\n", "                  Objective", "  Backward", "   Forward")
    printfmt("{1:19s} | {2:9s} {3:6s} | {4:6s} {5:6s} | {6:6s} {7:6s} |  Time\n", "      Expected", "  Bound", " % Gap", " Iters", " Time", " Iters", " Time")
end

#----------------------------------------------------------------------
# The following a modified version of that found at
#
# Humanize.jl    https://github.com/IainNZ/Humanize.jl
# Based on jmoiron's humanize Python library (MIT licensed):
#  https://github.com/jmoiron/humanize/
# All original code is (c) Iain Dunning and MIT licensed.
const gnu_suf = ["", "K", "M", "G", "T", "P", "E", "Z", "Y"]
function humanize(value::Int)
    if value < 1000 && value > -1000
        return humanize(value, "5d")
    else
        return humanize(value, "5.1f")
    end

end
function humanize(value::Number, format="5.1f")
    suffix  = gnu_suf
    base    = 1000.0
    bytes   = float(value)
    sig=sign(value)
    bytes   = abs(bytes)
    format  = "%$(format)%s"
    fmt_str = @eval (v,s)->@sprintf($format,v,s)
    unit    = base
    s       = suffix[1]
    for (i,s) in enumerate(suffix)
        unit = base ^ (i)
        bytes < unit && break
    end
    return fmt_str(sig*base * bytes / unit, s)
end
# End excerpt
#----------------------------------------------------------------------



function simulate(m::SDDPModel, n::Int, vars::Vector{Symbol}=[]; parallel=false)
    if parallel && length(workers()) < 2
        warn("Paralleisation requested but Julia is only running with a single processor. Start julia with `julia -p N` or use the `addprocs(N)` function to load N workers. Running in serial mode.")
        parallel = false
    end
    if parallel
        results = parallel_simulate(m, n, vars)
    else
        results = serial_simulate(m, n, vars)
    end
    results
end

end
