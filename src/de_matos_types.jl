"""
A single cut
"""
type Cut
    intercept::Float64
    coefficients::Vector{Float64}
end

"""
Cuts in a single stage problem
"""
type StageCuts
    n::Int                      # Number of cuts in stage problem
    X::Array{Float64, 2}        # x to evaluate at
    cuts::Vector{Cut}           # list of cuts in stage problem. length(cuts) == n
    nondominated::Vector{Int}   # number of points cut is nondominated at length(nondomindated) == n
    active_cut::Vector{Int}   # index of cut that is active at point x(i) length(active_cut) == size(X)[2]
end


function StageCuts(sp::Model, bound, n=1000)
    sd = stagedata(sp)

    X=rand(length(sd.state_vars),n)
    for (i,v) in enumerate(sd.state_vars)
        if getLower(v) == -Inf || getUpper(v) == Inf
            error("Cut deletion techniques can only be used if state variables contain non-infinite bounds. Currently, $(v) has a lower bound of $(getLower(v)) and an upper bound of $(getUpper(v)).")
        end
        X[i,:] = getLower(v) + (getUpper(v) - getLower(v)) * rand(n)
    end
    StageCuts(1, X, Cut[Cut(bound, zeros(size(X)[1]))], Int[size(X)[2]], ones(size(X)[2]))
end

"""
This function returns the active cut at x(i)
"""
getcut(stagecut::StageCuts, xi::Int) = stagecut.cuts[stagecut.active_cut[xi]]

function add_cut!(sense, cut::Cut, stagecut::StageCuts)
    stagecut.n += 1

    non_domination = size(stagecut.X)[2]
    for i=1:size(stagecut.X)[2]
        y0 = cut.intercept + dot(cut.coefficients, stagecut.X[:,i])
        c = getcut(stagecut, i)
        if is_dominated(sense, c.intercept + dot(c.coefficients, stagecut.X[:,i]), y0)
            non_domination -= 1
        else
            stagecut.nondominated[stagecut.active_cut[i]] -= 1
            stagecut.active_cut[i] = stagecut.n
        end
    end

    push!(stagecut.cuts, cut)
    push!(stagecut.nondominated, non_domination)

    return non_domination > 0
end

is_dominated(::Type{Val{:Min}}, y0::Float64, y1::Float64) = y0 < y1
is_dominated(::Type{Val{:Max}}, y0::Float64, y1::Float64) = y0 > y1
