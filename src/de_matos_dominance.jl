module deMatos

export Cut, ModelCuts

type Cut
    intercept::Float64
    coefficients::Vector{Float64}
end

type StageCuts
    n::Int
    X::Array{Float64, 2}
    cuts::Vector{Cut}
    nondominated::Vector{Int}
    active_cut::Vector{Int}
end
getcut(stagecuts::StageCuts, xi::Int) = stagecut.cuts[stagecut.active_cut[i]]

type ModelCuts
    cuts::Array{StageCuts, 2}
    level::Int
end

function add_cut!(cut::Cut, stagecut::StageCuts)
    stagecut.n += 1

    y0 = cut.intercept + dot(cut.coefficients, x)
    non_domination = size(X)[1]
    for i=1:size(X)[1]
        c = getcut(stagecut, i)
        if is_dominated(v, c.intercept + dot(c.coefficients, stagecut.X[:,i]), y0)
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

end

X = [[rand()] for i=1:10]

mc = deMatos.ModelCuts(X, 1, 1, 1)
deMatos.add_cut!(mc, deMatos.Cut(1., [1.]), 1, 1)
deMatos.add_cut!(mc, deMatos.Cut(1., [1.2]), 1, 1)
deMatos.add_cut!(mc, deMatos.Cut(2., [-1.]), 1, 1)

mc

println("The End")
