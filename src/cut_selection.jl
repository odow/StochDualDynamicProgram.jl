#  Copyright 2017, Oscar Dowson

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
    StageCuts(s.n, deepcopy(s.samplepoints), deepcopy(c.cuts), copy(s.nondominated), copy(s.activecut), s.writtentofile)
end

"""
    cutselection!(SDDPModel, CutSelectionMethod, iteration)

If `iteration` mod the CutSelectionMethod.frequency is zero, then the subproblems are rebuilt using the CutSelectionMethod.
"""
function cutselection!(m::SDDPModel, cutselection::CutSelectionMethod, iteration, print_level)
    if mod(iteration, cutselection.frequency) == 0
        print_level >= PRINTINFO && info("Running cut selection")
        rebuild_stageproblems!(m, cutselection)
    end
end
cutselection!(m::SDDPModel, cutselection::NoSelection, iteration) = nothing

# a wrapper for cut selection timings
function cutselection!(log::SolutionLog, m::SDDPModel, cutselection, iterations, print_level)
    if cutselection.frequency > 0
        tic()
        cutselection!(m, cutselection, iterations, print_level)
        log.time_cutselection += toq()
    end
    return
end

# true if y0 is dominated by y1
is_dominated(::Type{Min}, y0::Float64, y1::Float64) = y0 < y1
is_dominated(::Type{Max}, y0::Float64, y1::Float64) = y0 > y1

# add the cut 'cut' to the subproblem stagecut
function add_cut!{N}(sense::Sense, stagecut::StageCuts{N}, cut::Cut{N})
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

function rebuild_stageproblems!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, cutselection::CutSelectionMethod=NoSelection())
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
addsamplepoint!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, pass, t::Int, i::Int) = addsamplepoint!(X, stagecut(m, t, i), getx(m, pass, t))

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

function writecuts!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, cut_output_file::String)
    open(cut_output_file, "a") do f
        for t=1:(T-1)
            for i=1:M
                sc = stagecut(m, t, i)
                for cut in sc.cuts[(sc.writtentofile+1):end]
                    write(f, "$t, $i, $(cut.intercept)")
                    for coef in cut.coefficients
                        write(f, ", $(coef)")
                    end
                    write(f, "\n")
                end
                sc.writtentofile = length(sc.cuts)
            end
        end
    end
end

"""
    loadcuts!(model, filename)

This function loads cuts from the file at `filename` to the SDDPModel `model`.
"""
function loadcuts!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, filename::String)
    open(filename, "r") do f
        while true
            line = readline(f)
            (line == nothing || line == "") && break
            line = split(strip(line), ",")
            @assert length(line) >= 3
            t = parse(Int, line[1])
            i = parse(Int, line[2])
            theta = parse(Float64, line[3])
            if length(line) > 4
                xcoeff = map(x->parse(Float64, x), line[4:end])
            else
                xcoeff = [parse(Float64, line[4])]
            end
            cut = Cut(theta, xcoeff)
            # add to cutselection storage
            add_cut!(X, stagecut(m, t, i), cut)
            # add to problem
            addcut!(X, subproblem(m, t, i), cut)
        end
    end
end

# This function recalculates dominance
function recalculate_dominance!{N}(sense::Sense, sc::StageCuts{N})
    sc.n = length(sc.samplepoints)
    sc.activecut = zeros(sc.n)
    sc.nondominated = zeros(length(sc.cuts))
    for (sample_idx, samplepoint) in enumerate(sc.samplepoints)
        best_idx = 1
        best_y = evaluate(sc.cuts[1], samplepoint)
        for cut_idx in 2:length(sc.cuts)
            y = evaluate(sc.cuts[cut_idx], samplepoint)
            if is_dominated(sense, best_y, y)
                best_idx = cut_idx
                best_y = y
            end
        end
        sc.nondominated[best_idx] += 1
        sc.activecut[sample_idx] = best_idx
    end
end
recalculate_dominance!(::CutSelectionMethod, m) = nothing
function recalculate_dominance!{T, M, S, X, TM}(::LevelOne, m::SDDPModel{T, M, S, X, TM})
    for t=1:T
        for i=1:M
            recalculate_dominance!(X, stagecut(m, t, i))
        end
    end
end

# Run the exact method for cut selection
function rebuildcuts!{N}(::Deterministic, sense::Sense, sp::Model, sc::StageCuts{N}) #, bound::Real=-Inf)
    m = Model()
    x = @variable(m, [i=1:N], lowerbound=getlowerbound(stagedata(sp).state_vars[i]), upperbound=getupperbound(stagedata(sp).state_vars[i]))
    y = @variable(m)
    for cut in sc.cuts
        if isa(sense, Max)
            @constraints(m, begin
                # y <= bound
                y <= dot(cut, x)
            end)
        else
            @constraints(m, begin
                # y >= bound
                y >= dot(cut, x)
            end)
        end
    end
    activecuts = Cut{N}[]
    for cut in sc.cuts
        if isa(sense, Max)
            @objective(m, Min, dot(cut, x) - y)
        else
            @objective(m, Min, y - dot(cut, x))
        end
        solve(m)
        if getobjectivevalue(m) < 0.
            push!(activecuts, cut)
        end
    end
    sc.activecut = activecuts
    recalculate_dominance!(sense, sc)
end

Base.dot(cut::Cut, x::Vector{JuMP.Variable}) = cut.intercept + dot(cut.coefficients, x)
