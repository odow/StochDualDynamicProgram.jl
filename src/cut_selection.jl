# Overload base.dot to handle (vector x NTuple)
function Base.dot{T<:Real, N}(x::Vector{T}, y::NTuple{N, T})
    @assert length(x) == N
    z = zero(T)
    @inbounds for i=1:N
        z += x[i] * y[i]
    end
    z
end
Base.dot{T<:Real, N}(x::NTuple{N, T}, y::Vector{T}) = dot(y,x)

# Define some function to allow comparison of cuts
Base.hash{N}(c::Cut{N}) = hash(c.intercept, hash(c.coefficients))
Base.isequal{N}(c1::Cut{N}, c2::Cut{N}) = c1.intercept == c2.intercept && c1.coefficients == c2.coefficients
Base.copy{N}(c::Cut{N}) = Cut{N}(c.intercept, copy(c.coefficients))

function Base.copy{N}(s::StageCuts{N})
    StageCuts(s.n, deepcopy(s.samplepoints), deepcopy(c.cuts), copy(s.nondominated), copy(s.activecut))
end

# true if y0 is dominated by y1
is_dominated(::Type{Min}, y0::Float64, y1::Float64) = y0 < y1
is_dominated(::Type{Max}, y0::Float64, y1::Float64) = y0 > y1

# add the cut 'cut' to the subproblem stagecut
function add_cut!{N}(sense::Sense, cut::Cut{N}, stagecut::StageCuts{N})
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

# Evalute the cut c at x position t
evaluate{N}(c::Cut{N}, t::NTuple{N, Float64}) = c.intercept + dot(c.coefficients, t)

# This function returns the active cut at x(i)
function getactivecut(stagecut::StageCuts, xi::Int)
    @assert xi <= stagecut.n
    stagecut.cuts[stagecut.activecut[xi]]
end

function rebuild_stageproblems!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, cutselection::CutSelectionMethod)
    create_subproblems!(m)
    cuts_added = 0
    total_cuts = 0
    for t=1:(T-1)
        for i=1:M
            sp = subproblem(m, t, i)
            cuts_added += rebuildcuts!(cutselection, X, sp, stagecut(m, t, i))
            total_cuts += length(stagecut(m, t, i).cuts)
        end
    end
end

function addsamplepoint!{N}(sense::Sense, stagecut::StageCuts{N}, x::Vector{Float64})
    @assert length(x) == N
    # add sample point
    tup = tuple(copy(x)...)
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
addsamplepoint!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t::Int, i::Int) = addsamplepoint!(X, stagecut(m, t, i), getx(m, t, i))

function rebuildcuts!(::NoSelection, sense, sp, stagecut)
    cuts_added = 0
    for cut in unique(stagecut.cuts)
        addcut!(sense, sp, cut)
        cuts_added += 1
    end
    cuts_added
end

function rebuildcuts!(::LevelOne, sense, sp, stagecut)
    cuts_added = 0
    for xi in unique(stagecut.activecut)
        addcut!(sense, sp, stagecut.cuts[xi])
        cuts_added += 1
    end
    cuts_added
end

# # This function recalculates dominance
# function recalculate_dominance!{N}(sense::Sense, sc::StageCuts{N})
#     sc.n = length(sc.samplepoints)
#     sc.activecut = zeros(sc.n)
#     sc.nondominated = zeros(length(sc.cuts))
#     for (sample_idx, samplepoint) in enumerate(sc.samplepoints)
#         best_idx = 1
#         best_y = evaluate(sc.cuts[1], samplepoint)
#         for cut_idx in 2:length(sc.cuts)
#             y = evaluate(sc.cuts[cut_idx], samplepoint)
#             if is_dominated(sense, best_y, y)
#                 best_idx = cut_idx
#                 best_y = y
#             end
#         end
#         sc.nondominated[best_idx] += 1
#         sc.activecut[sample_idx] = best_idx
#     end
# end
#
#
# function rebuild_stageproblems!(::Deterministic, m::SDDPModel)
#     deterministic_prune!(m)
#     rebuild_stageproblems!(NoSelection(), m)
# end
# rebuild_stageproblems!(m::SDDPModel) = rebuild_stageproblems!(NoSelection(), m)
#
# function add_cut!(sense, sp::Model, cut::Cut, method::CutSelectionMethod=NoSelection())
#     @expression(sp, rhs, cut.intercept + sum{coeff * stagedata(sp).state_vars[i], (i, coeff) in enumerate(cut.coefficients)})
#     add_cut!(sense, sp, rhs, method)
# end
#
# # Run the exact method for cut selection
# function deterministic_prune!{N}(sense::Sense, bound::Real, sp::Model, sc::StageCuts{N})
#     m = Model()
#     @variable(m, getlowerbound(stagedata(sp).state_vars[i]) <= x[i=1:N] <= getupperbound(stagedata(sp).state_vars[i]))
#     @variable(m, y)
#     for cut in sc.cuts
#         if isa(sense, Max)
#             @constraints(m, begin
#                 y <= bound
#                 y <= cut.intercept + sum{cut.coefficients[i] * x[i], i=1:N}
#             end)
#         else
#             @constraints(m, begin
#                 y >= bound
#                 y >= cut.intercept + sum{cut.coefficients[i] * x[i], i=1:N}
#             end)
#         end
#     end
#     activecuts = Cut{N}[]
#     for cut in sc.cuts
#         if isa(sense, Max)
#             @objective(m, Min, cut.intercept + sum{cut.coefficients[i] * x[i], i=1:N} - y)
#         else
#             @objective(m, Min, y - cut.intercept + sum{cut.coefficients[i] * x[i], i=1:N})
#         end
#         solve(m)
#         if getobjectivevalue(m) < 0.
#             push!(activecuts, cut)
#         end
#     end
#     sc.activecut = activecuts
#     recalculate_dominance!(sense, sc)
# end
#
# # This function runs the deterministic prune on the entire SDDPModel
# function deterministic_prune!(m::SDDPModel)
#     for i=1:m.stages
#         for j=1:m.markov_states
#             deterministic_prune!(m.sense, m.value_to_go_bound, m.stage_problems[i,j], m.stagecuts[i,j])
#         end
#     end
# end
#
# """
# This function loads cuts from a file.
# """
# function loadcuts!(m::SDDPModel, filename::ASCIIString)
#     open(filename, "r") do f
#         while true
#             line = readline(f)
#             (line == nothing || line == "") && break
#             line = split(strip(line), ",")
#             @assert length(line) >= 3
#             stage = parse(Int, line[1])
#             markov_state = parse(Int, line[2])
#             sp = m.stage_problems[stage, markov_state]
#             theta = parse(Float64, line[3])
#
#             @assert length(line) == (3 + length(stagedata(sp).state_vars))
#             if length(line) > 4
#                 xcoeff = map(x->parse(Float64, x), line[4:end])
#             else
#                 xcoeff = [parse(Float64, line[4])]
#             end
#             loadcut!(m.sense, sp, theta, xcoeff)
#         end
#     end
# end
# function loadcut!(::Type{Max}, sp::Model, theta::Float64, xcoeff::Vector{Float64})
#     @constraint(sp, stagedata(sp).theta <= theta + sum{xcoeff[i] * v, (i, v) in enumerate(stagedata(sp).state_vars)})
# end
# function loadcut!(::Type{Min}, sp::Model, theta::Float64, xcoeff::Vector{Float64})
#     @constraint(sp, stagedata(sp).theta >= theta + sum{xcoeff[i] * v, (i, v) in enumerate(stagedata(sp).state_vars)})
# end
