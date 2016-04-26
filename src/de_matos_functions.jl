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
function rebuild_stageproblems!(m::SDDPModel)
    create_subproblems!(m)
    cuts_added = 0
    total_cuts = 0
    for stage=1:(m.stages-1)
        for markovstate=1:m.markov_states
            sp = m.stage_problems[stage, markovstate]
            stagecut = m.stagecuts[stage, markovstate]

            for xi in unique(stagecut.activecut)
                add_cut!(m.sense, sp, stagecut.cuts[xi])
                cuts_added += 1
            end
            total_cuts += length(stagecut.cuts)
        end
    end
    # info("Rebuilding model using $(cuts_added) of $(total_cuts) ($(round(cuts_added/total_cuts*100, 2))\%) discovered cuts.")
end

function add_cut!(sense, sp::Model, cut::Cut)
    @defExpr(rhs, cut.intercept + sum{coeff * stagedata(sp).state_vars[i], (i, coeff) in enumerate(cut.coefficients)})
    add_cut!(sense, sp, rhs)
end
