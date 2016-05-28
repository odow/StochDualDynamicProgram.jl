function textify(s::SolutionLog)
    string(s.ci_lower, ",", s.ci_upper, ",", s.bound, ",", s.cuts, ",", s.time_backwards, ",", s.simulations, ",", s.time_forwards, ",", s.time_cutselection, "\n")
end

function log!(s::Solution, ci_lower::Float64, ci_upper::Float64, bound::Float64, cuts::Int, time_backwards::Float64, simulations::Int, time_forwards::Float64, time_cutselection::Float64)
    push!(s.trace, SolutionLog(ci_lower, ci_upper, bound, cuts, time_backwards, simulations, time_forwards, time_cutselection))
end

function getTrace(sol::Solution, sym::Symbol)
    if !(sym in fieldnames(sol))
        error("[sym] in getTrace(sol::Solution, sym::Symbol) must be one of $(fieldnames(sol))")
    end
    return [s.(sym) for s in sol.trace]
end

function Base.print(m::SDDPModel, l::SolutionLog, filename::ASCIIString)
    open(filename, "a") do f
        write(f, textify(l))
    end
    print(m, l)
end
Base.print(m::SDDPModel, l::SolutionLog, filename::Void) = print(m, l)
function Base.print(m::SDDPModel, l::SolutionLog)
    printfmt("{1:>9s} {2:>9s} | {3:>9s} {4:>6.2f} | {5:6s} {6:6s} | {7:6s} {8:6s} | {9:5s}\n",
        humanize(l.ci_lower, "7.3f"), humanize(l.ci_upper, "7.3f"), humanize(l.bound, "7.3f"), 100*rtol(m), humanize(l.cuts), humanize(l.time_backwards), humanize(l.simulations), humanize(l.time_forwards), humanize(l.time_cutselection, "6.2f"))
end

getCloseCIBound(::Type{Min}, m::SDDPModel) = m.confidence_interval[1]
getCloseCIBound(::Type{Max}, m::SDDPModel) = m.confidence_interval[2]
getCloseCIBound{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}) = getCloseCIBound(X, m)
getBound(m::SDDPModel) = m.valid_bound
atol(::Type{Min}, m::SDDPModel) = getCloseCIBound(m) - getBound(m)
atol(::Type{Max}, m::SDDPModel) = getBound(m) - getCloseCIBound(m)
atol{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}) = atol(X, m)

"""
    rtol(SDDPModel)

Relative tolerance of the solution
Defined as [Outer bound - closest simulated bound] / [Outer bound]
"""
function rtol(m)
    abs(getBound(m)) != Inf?atol(m) / abs(getBound(m)):Inf
end

function printheader()
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
