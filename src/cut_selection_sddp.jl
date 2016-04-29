function initialise_cut_selection!(m::SDDPModel)
    m.stagecuts = Array(StageCuts, (m.stages, m.markov_states))
    for stage=1:m.stages
        for markov_state=1:m.markov_states
            m.stagecuts[stage,markov_state] = StageCuts(m.stage_problems[stage,markov_state], m.value_to_go_bound)
        end
    end
end

function add_cut!(sense, sp::Model, rhs::JuMP.GenericAffExpr, stagecut::StageCuts)
    add_cut!(sense, Cut(rhs.constant, aggregate_terms(sp, rhs)), stagecut)
end

"""
This function rebuilds the stageproblems using the nondominated cuts
"""
function add_cuts!(::NoSelection, sense, sp, stagecut)
    cuts_added = 0
    for cut in unique(stagecut.cuts)
        add_cut!(sense, sp, cut)
        cuts_added += 1
    end
    cuts_added
end

function add_cuts!(::LazyConstraint, sense, sp, stagecut)
    cuts_added = 0
    for cut in unique(stagecut.cuts)
        add_cut!(sense, sp, cut, LazyConstraint())
        cuts_added += 1
    end
    cuts_added
end

function add_cuts!(::LevelOne, sense, sp, stagecut)
    cuts_added = 0
    for xi in unique(stagecut.activecut)
        add_cut!(sense, sp, stagecut.cuts[xi])
        cuts_added += 1
    end
    cuts_added
end

function rebuild_stageproblems!(valswitch, m::SDDPModel)
    create_subproblems!(m)
    cuts_added = 0
    total_cuts = 0
    for stage=1:(m.stages-1)
        for markovstate=1:m.markov_states
            sp = m.stage_problems[stage, markovstate]
            stagecut = m.stagecuts[stage, markovstate]
            cuts_added += add_cuts!(valswitch, m.sense, sp, stagecut)
            total_cuts += length(stagecut.cuts)
        end
    end
    # info("Rebuilding model using $(cuts_added) of $(total_cuts) ($(round(cuts_added/total_cuts*100, 2))\%) discovered cuts.")
end

function rebuild_stageproblems!(::Deterministic, m::SDDPModel)
    deterministic_prune!(m)
    rebuild_stageproblems!(NoSelection(), m)
end
rebuild_stageproblems!(m::SDDPModel) = rebuild_stageproblems!(NoSelection(), m)

function add_cut!(sense, sp::Model, cut::Cut, method::CutSelectionMethod=NoSelection())
    @defExpr(rhs, cut.intercept + sum{coeff * stagedata(sp).state_vars[i], (i, coeff) in enumerate(cut.coefficients)})
    add_cut!(sense, sp, rhs, method)
end

function deterministic_prune!(m::SDDPModel)
    for i=1:m.stages
        for j=1:m.markov_states
            deterministic_prune!(m.sense, m.value_to_go_bound, m.stage_problems[i,j], m.stagecuts[i,j])
        end
    end
end

function deterministic_prune!{N}(sense, bound, sp::Model, sc::StageCuts{N})
    m = Model()
    @defVar(m, getLower(stagedata(sp).state_vars[i]) <= x[i=1:N] <= getUpper(stagedata(sp).state_vars[i]))
    @defVar(m, y)
    for cut in sc.cuts
        if isa(sense, Val{:Max})
            @addConstraint(m, y <= bound)
            @addConstraint(m, y <= cut.intercept + sum{cut.coefficients[i] * x[i], i=1:N})
        else
            @addConstraint(m, y >= bound)
            @addConstraint(m, y >= cut.intercept + sum{cut.coefficients[i] * x[i], i=1:N})
        end
    end
    activecuts = Cut{N}[]
    for cut in sc.cuts
        if isa(sense, Val{:Max})
            @setObjective(m, Min, cut.intercept + sum{cut.coefficients[i] * x[i], i=1:N} - y)
        else
            @setObjective(m, Min, y - cut.intercept + sum{cut.coefficients[i] * x[i], i=1:N})
        end
        solve(m)
        if getObjectiveValue(m) < 0.
            push!(activecuts, cut)
        end
    end
    sc.activecut = activecuts
    recalculate_dominance!(sense, sc)
end
