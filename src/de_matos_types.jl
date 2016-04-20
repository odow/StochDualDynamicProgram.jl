import Base.dot
function Base.dot{T<:Real, N}(x::Vector{T}, y::NTuple{N, T})
    @assert length(x) == N
    z = zero(T)
    @inbounds for i=1:N
        z += x[i] * y[i]
    end
    z
end
Base.dot{T<:Real, N}(x::NTuple{N, T}, y::Vector{T}) = dot(y,x)

"""
A single cut
"""
type Cut
    intercept::Float64
    coefficients::Vector{Float64}
    Cut(intercept, coefficients::Vector) = new(round(intercept, 10), round(coefficients, 10))
end


"""
Cuts in a single stage problem
"""
type StageCuts{N}
    n::Int                      # Number of sample points in stage problem
    samplepoints::Vector{NTuple{N, Float64}}  # x to evaluate at
    cuts::Vector{Cut}           # list of cuts in stage problem.
    nondominated::Vector{Int}   # number of points cut is nondominated at length(nondomindated) == length(cuts)
    activecut::Vector{Int}      # index of cut that is active at point x(i) length(active_cut) == n
end
StageCuts(sp::Model, bound) = StageCuts(0, NTuple{length(stagedata(sp).state_vars), Float64}[], Cut[Cut(bound, zeros(length(stagedata(sp).state_vars)))], Int[0], Int[])

function addsamplepoint!(sense, stagecut::StageCuts, x::Vector{Float64})
    # add sample point
    tup = tuple(round(copy(x), 10)...)
    if !(tup in stagecut.samplepoints)
        stagecut.n += 1
        push!(stagecut.samplepoints, tup)

        # calculate nondominated cut
        y0 = evaluate(stagecut.cuts[1], tup)
        iBest = 1
        for i=2:length(stagecut.cuts)
            y1 = evaluate(stagecut.cuts[i], tup)
            if is_dominated(sense, y0, y1)
                y0 = y1
                iBest = i
            end
        end
        push!(stagecut.activecut, iBest)
        stagecut.nondominated[iBest] += 1
    end

    return
end

evaluate{N}(c::Cut, t::NTuple{N, Float64}) = c.intercept + dot(c.coefficients, t)

"""
This function returns the active cut at x(i)
"""
function getcut(stagecut::StageCuts, xi::Int)
    @assert xi <= stagecut.n
    stagecut.cuts[stagecut.activecut[xi]]
end

function add_cut!(sense, cut::Cut, stagecut::StageCuts)
    push!(stagecut.cuts, cut)

    non_domination = stagecut.n
    for i=1:stagecut.n
        y0 = evaluate(cut, stagecut.samplepoints[i])
        c = getcut(stagecut, i)
        if is_dominated(sense, evaluate(c, stagecut.samplepoints[i]), y0)
            non_domination -= 1
        else
            stagecut.nondominated[stagecut.activecut[i]] -= 1
            stagecut.activecut[i] = length(stagecut.cuts)
        end
    end

    push!(stagecut.nondominated, non_domination)

    return non_domination > 0
end

is_dominated(::Type{Val{:Min}}, y0::Float64, y1::Float64) = y0 <= y1
is_dominated(::Type{Val{:Max}}, y0::Float64, y1::Float64) = y0 >= y1
