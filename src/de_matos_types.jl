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
type Cut{N}
    intercept::Float64
    coefficients::Vector{Float64}
end
Cut(intercept, coefficients::Vector) = Cut{length(coefficients)}(round(intercept, 15), round(coefficients, 15))
Base.hash{N}(c::Cut{N}) = hash(c.intercept, hash(c.coefficients))
Base.isequal{N}(c1::Cut{N}, c2::Cut{N}) = c1.intercept == c2.intercept && c1.coefficients == c2.coefficients
Base.copy{N}(c::Cut{N}) = Cut{N}(c.intercept, copy(c.coefficients))

"""
Cuts in a single stage problem
"""
type StageCuts{N}
    n::Int                      # Number of sample points in stage problem
    samplepoints::Vector{NTuple{N, Float64}}  # x to evaluate at
    cuts::Vector{Cut{N}}           # list of cuts in stage problem.
    nondominated::Vector{Int}   # number of points cut is nondominated at length(nondomindated) == length(cuts)
    activecut::Vector{Int}      # index of cut that is active at point x(i) length(active_cut) == n
end
function StageCuts(sp::Model, bound)
    N = length(stagedata(sp).state_vars)
    StageCuts(0, NTuple{N, Float64}[], Cut{N}[Cut(bound, zeros(N))], Int[0], Int[])
end
function Base.copy{N}(s::StageCuts{N})
    StageCuts(s.n, deepcopy(s.samplepoints), deepcopy(c.cuts), copy(s.nondominated), copy(s.activecut))
end

function addsamplepoint!{N}(sense, stagecut::StageCuts{N}, x::Vector{Float64})
    @assert length(x) == N

    # add sample point
    tup = tuple(round(copy(x), 15)...)

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

evaluate{N}(c::Cut{N}, t::NTuple{N, Float64}) = c.intercept + dot(c.coefficients, t)

"""
This function returns the active cut at x(i)
"""
function getactivecut(stagecut::StageCuts, xi::Int)
    @assert xi <= stagecut.n
    stagecut.cuts[stagecut.activecut[xi]]
end

function add_cut!{N}(sense, cut::Cut{N}, stagecut::StageCuts{N})
    push!(stagecut.cuts, cut)

    non_domination = stagecut.n
    for i=1:stagecut.n
        y_new = evaluate(cut, stagecut.samplepoints[i])
        c = getactivecut(stagecut, i)
        if is_dominated(sense, y_new, evaluate(c, stagecut.samplepoints[i]))
            non_domination -= 1
        else
            stagecut.nondominated[stagecut.activecut[i]] -= 1
            stagecut.activecut[i] = length(stagecut.cuts)
        end
    end

    push!(stagecut.nondominated, non_domination)

    return non_domination > 0
end

# true if y0 is dominated by y1
is_dominated(::Type{Val{:Min}}, y0::Float64, y1::Float64) = y0 < y1
is_dominated(::Type{Val{:Max}}, y0::Float64, y1::Float64) = y0 > y1
