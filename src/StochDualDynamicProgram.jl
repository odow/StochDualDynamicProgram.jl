module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp
using Formatting
using Distributions, StatsBase

export SDDPModel,
    @defStateVar, @defValueToGo, @addScenarioConstraint, @setStageProfit,
    simulate, load_cuts!,
    LevelOne, Deterministic, NoSelection,
    ForwardPass, BackwardPass, Parallel,
    NestedCVar,
    Convergence

const CONVERGENCE_TERMINATION = :Convergenced
const ITERATION_TERMINATION = :MaximumIterations
const UNKNOWN_TERMINATION = :Unknown

abstract RiskMeasure
immutable NestedCVar <: RiskMeasure
    beta::Float64
    lambda::Float64
    NestedCVar{T<:Real, S<:Real}(beta::T, lambda::S) = (@assert beta >= 0. && beta <= 1. && lambda >= 0. && lambda <= 1.; new(beta, lambda))
end
NestedCVar(;beta=1., lambda=1.) = NestedCVar(beta, lambda)
Expectation() = NestedCVar(1., 1.)

type Convergence
    n
    frequency::Int
    terminate::Bool
    quantile::Float64
    function Convergence(simulations, frequency::Int, terminate::Bool=false, quantile::Float64=0.95)
        @assert quantile >= 0 && quantile <= 1.
        new(simulations, frequency, terminate, quantile)
    end
end
Convergence(;simulations=1, frequency=1, terminate=false, quantile=0.95) = Convergence(simulations, frequency, terminate, quantile)

include("macros.jl")
include("cut_selection.jl")
include("SDDPModel.jl")
include("cut_selection_sddp.jl")
include("SDDPalgorithm.jl")
include("parallel.jl")

type SolutionLog
    ci_lower::Float64
    ci_upper::Float64
    bound::Float64
    cuts::Int
    time_backwards::Float64
    time_forwards::Float64
end

type Solution
    status::Symbol
    trace::Vector{SolutionLog}
end
Solution() = Solution(UNKNOWN_TERMINATION, SolutionLog[])
setStatus!(s::Solution, sym::Symbol) = (s.status = sym)

function log!(s::Solution, ci_lower::Float64, ci_upper::Float64, bound::Float64, cuts::Int, time_backwards::Float64, time_forwards::Float64)
    push!(s.trace, SolutionLog(ci_lower, ci_upper, bound, cuts, time_backwards, time_forwards))
end

function getTrace(sol::Solution, sym::Symbol)
    if !(sym in fieldnames(sol))
        error("[sym] in getTrace(sol::Solution, sym::Symbol) must be one of $(fieldnames(sol))")
    end
    return [s.(sym) for s in sol.trace]
end

"""
Solve the model using the SDDP algorithm.
"""
function JuMP.solve(m::SDDPModel; maximum_iterations=1, convergence=Convergence(1, 1, false, 0.95), risk_measure=Expectation(), cut_selection=NoSelection(), parallel=Parallel(), simulation_passes=-1, log_frequency=-1)

    if simulation_passes != -1 || log_frequency != -1
        warn("The options [log_frequency] and [simulation_passes] are deprecated. Use [convergence=Convergence(simulation_passes, frequency::Int, terminate::Bool)] instead.")
        log_frequency != -1 && (convergence.frequency = log_frequency)
        simulation_passes != -1 && (convergence.n = simulation_passes)
    end

    @assert maximum_iterations >= 0

    if !isa(parallel, Parallel{Serial, Serial}) && length(workers()) < 2
        warn("Paralleisation requested (cuts_per_processor=$(cuts_per_processor)) but Julia is only running with a single processor. Start julia with `julia -p N` or use the `addprocs(N)` function to load N workers. Running in serial mode.")
        parallel = Parallel()
    end

    # Initial output for user
    print_stats_header()

    solution = Solution()

    m.QUANTILE = convergence.quantile

    # try
    solve!(m, solution, convergence, maximum_iterations, risk_measure, cut_selection, parallel)
    # catch InterruptException
        # warn("Terminating early")
    # end
    solution
end

terminate(converged::Bool, do_test::Bool) = converged & do_test
terminate(converged::Bool, do_test::Bool, simulation_passes::Range) = converged
terminate(converged::Bool, do_test::Bool, simulation_passes) = terminate(converged, do_test)

forward_pass!(ty::Serial, m::SDDPModel, simulation_passes) = forward_pass!(m, simulation_passes)
forward_pass!(ty::ForwardPass, m::SDDPModel, simulation_passes) = parallel_forward_pass!(m, simulation_passes)

backward_pass!(ty::BackwardPass, m::SDDPModel, method::CutSelectionMethod) = parallel_backward_pass!(m, ty.cuts_per_processor, method)
backward_pass!(ty::Serial, m::SDDPModel, method::CutSelectionMethod) = backward_pass!(m, method.frequency > 0, method)

function solve!(m::SDDPModel, solution::Solution, convergence::Convergence, maximum_iterations::Int, risk_measure::RiskMeasure, cut_selection::CutSelectionMethod, parallel::Parallel)

    # Intialise model on worker processors
    if !isa(parallel, Parallel{Serial, Serial})
        nworkers = initialise_workers!(m)
    end

    # Initialise cut selection storage since we need it to communicate between processors
    initialise_cut_selection!(m)

    # Set risk aversion parameters
    m.beta_quantile, m.risk_lambda = risk_measure.beta, risk_measure.lambda

    iterations=0
    nsimulations = 0
    last_convergence_test, last_cutselection = 0, 0
    time_backwards, time_forwards, time_cutselection = 0., 0., 0.
    npass = isa(parallel.backward_pass, BackwardPass)?(parallel.backward_pass.cuts_per_processor * nworkers):1

    while iterations < maximum_iterations

        # Update iteration counters
        iterations += npass
        last_convergence_test += npass
        last_cutselection += npass

        # Cutting passes
        tic()
        backward_pass!(parallel.backward_pass, m, cut_selection)
        time_backwards += toq()

        # Rebuild models if using Cut Selection
        if cut_selection.frequency > 0 && last_cutselection >= cut_selection.frequency
            tic()
            rebuild_stageproblems!(cut_selection, m)
            time_cutselection += toq()
        end

        # Test convergence if appropriate
        if convergence.frequency > 0 && last_convergence_test >= convergence.frequency
            last_convergence_test = 0
            tic()
            (is_converged, n) = forward_pass!(parallel.forward_pass, m, convergence.n)
            nsimulations += Int(n)
            time_forwards += toq()

            # Output to user
            log!(solution, m.confidence_interval[1], m.confidence_interval[2], m.valid_bound, iterations, time_backwards, time_forwards)
            print_stats(m, iterations, time_backwards, nsimulations, time_forwards, time_cutselection)

            # Terminate if converged and appropriate
            if terminate(is_converged, convergence.terminate, convergence.n)
                # We have terminated due to convergence test
                setStatus!(solution, CONVERGENCE_TERMINATION)
                return
            end
        end
    end
    # We have terminated due to iteration limit
    setStatus!(solution, ITERATION_TERMINATION)
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
