function add_cut!(sense, sp::Model, rhs::JuMP.GenericAffExpr, stagecut::StageCuts)
    add_cut!(sense, Cut(rhs.constant, aggregate_terms(sp, rhs)), stagecut)
end

"""
This function rebuilds the stageproblems using the nondominated cuts
"""
function rebuild_stageproblems!(m::SDDPModel)
    create_subproblems!(m)
    for stage=1:(m.stages-1)
        for markovstate=1:m.markov_states
            sp = m.stage_problems[stage, markovstate]
            stagecut = m.stagecuts[stage, markovstate]
            for xi in unique(stagecut.activecut)
                add_cut!(m.sense, sp, getcut(stagecut, xi))
            end
        end
    end
end

function add_cut!(sense, sp::Model, cut::Cut)
    @defExpr(rhs, cut.intercept + sum{coeff * stagedata(sp).state_vars[i], (i, coeff) in enumerate(cut.coefficients)})
    add_cut!(sense, sp, rhs)
end
