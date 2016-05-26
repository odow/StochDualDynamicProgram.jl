const CONVERGENCE_TERMINATION = :Convergenced
const ITERATION_TERMINATION = :MaximumIterations
const UNKNOWN_TERMINATION = :Unknown

# ==============================================================================
#   Optimisation Sense
abstract Sense
immutable Max <: Sense end
immutable Min <: Sense end

# ==============================================================================
#   Risk Measures
abstract RiskMeasure

immutable NestedCVar <: RiskMeasure
    beta::Float64
    lambda::Float64
    NestedCVar{T<:Real, S<:Real}(beta::T, lambda::S) = (@assert beta >= 0. && beta <= 1. && lambda >= 0. && lambda <= 1.; new(beta, lambda))
end

NestedCVar(;beta=1., lambda=1.) = NestedCVar(beta, lambda)
Expectation() = NestedCVar(1., 1.)

# ==============================================================================
#   Convergence
type Convergence
    n
    frequency::Int
    terminate::Bool
    quantile::Float64
    variancereduction::Bool
    function Convergence(simulations, frequency::Int, terminate::Bool=false, quantile::Float64=0.95, variancereduction::Bool=true)
        @assert quantile >= 0 && quantile <= 1.
        new(simulations, frequency, terminate, quantile, variancereduction)
    end
end

Convergence(;simulations=1, frequency=1, terminate=false, quantile=0.95, variancereduction=true) = Convergence(simulations, frequency, terminate, quantile, variancereduction)

# ==============================================================================
#   Regularisation
abstract Regularisation

immutable NoRegularisation <: Regularisation
    initial::Float64
    decayrate::Float64
    NoRegularisation() = new(0., 0.)
end

immutable LinearRegularisation <: Regularisation
    initial::Float64
    decayrate::Float64
    LinearRegularisation(initial=1., decayrate=0.95) = new(initial, decayrate)
end

immutable QuadraticRegularisation <: Regularisation
    initial::Float64
    decayrate::Float64
    QuadraticRegularisation(initial=1., decayrate=0.95) = new(initial, decayrate)
end

# ==============================================================================
#   Cut Selection
abstract CutSelectionMethod

immutable LevelOne <: CutSelectionMethod
    frequency::Int
    LevelOne(frequency::Int) = (@assert frequency > 0; new(frequency))
end

immutable Deterministic <: CutSelectionMethod
    frequency::Int
    Deterministic(frequency::Int) = (@assert frequency > 0; new(frequency))
end

immutable NoSelection <: CutSelectionMethod
    frequency::Int
    NoSelection() = new(0)
end

# ==============================================================================
#   Solution
type SolutionLog
    ci_lower::Float64
    ci_upper::Float64
    bound::Float64
    cuts::Int
    time_backwards::Float64
    simulations::Int
    time_forwards::Float64
    time_cutselection::Float64
end
function textify(s::SolutionLog)
    string(s.ci_lower, ",", s.ci_upper, ",", s.bound, ",", s.cuts, ",", s.time_backwards, ",", s.simulations, ",", s.time_forwards, ",", s.time_cutselection, "\n")
end

type Solution
    status::Symbol
    trace::Vector{SolutionLog}
end
Solution() = Solution(UNKNOWN_TERMINATION, SolutionLog[])

setStatus!(s::Solution, sym::Symbol) = (s.status = sym)

function log!(s::Solution, ci_lower::Float64, ci_upper::Float64, bound::Float64, cuts::Int, time_backwards::Float64, simulations::Int, time_forwards::Float64, time_cutselection::Float64)
    push!(s.trace, SolutionLog(ci_lower, ci_upper, bound, cuts, time_backwards, simulations, time_forwards, time_cutselection))
end

function getTrace(sol::Solution, sym::Symbol)
    if !(sym in fieldnames(sol))
        error("[sym] in getTrace(sol::Solution, sym::Symbol) must be one of $(fieldnames(sol))")
    end
    return [s.(sym) for s in sol.trace]
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
