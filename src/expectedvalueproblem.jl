# We solve the expected value problem first to obtain a set of cuts that roughly
# define the cost function. We hope that this improves convergence.
Base.dot(x::WeightVec, y::Vector) = dot(x.values, y)
Base.dot(x::Vector, y::WeightVec) = dot(y, x)

function setexpectedvaluerhs!(m::SDDPModel, sp::Model)
    # @assert is_subproblem(sp)
    for (constr, Ω) in stagedata(sp).scenario_constraints
        JuMP.setRHS(constr, dot(m.scenario_probability, Ω))
    end
end

function expectedforwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int, cutselection::CutSelectionMethod, forwardpass::ForwardPass)
    markov = 0                                       # initialise
    for pass = 1:n                                   # for n passes
        markov = m.initial_markov_state              # initial markov state
        if m.initial_markov_state==0                 # no initial state specified
            markov = transition(m, 1, markov)        # transition immediately
        end
        m.forwardstorage.obj[pass] = 0               # reset objective
        for t=1:T                                    # for each stage
            savemarkov!(m, pass, t, markov)          # store markov state
            sp = subproblem(m, t, markov)            # get subproblem

            setexpectedvaluerhs!(m, sp)
            regularisedsolve!(X, sp, forwardpass.regularisation)

            saveobj!(m, pass, getstagevalue(sp))     # store objective
            savex!(m, sp, pass, t)                   # store state
            addsamplepoint!(m, cutselection,         # store sample points for
                    pass, t, markov)                 #    cutselection
            if t < T                                 # don't do this for the last stage
                pass_states!(m, sp, t, forwardpass.regularisation) # pass state values forward
                markov = transition(m, t, markov, forwardpass.uniformsampling)    # transition
            end
        end
    end
end

function expectedforwardpass!{T, M, S, X, TM}(log, m::SDDPModel{T, M, S, X, TM}, n::Int, cutselection::CutSelectionMethod, forwardpass::ForwardPass)
    tic()
    expectedforwardpass!(m, n, cutselection, forwardpass)
    log.simulations += n
    log.time_forwards += toq()
    return
end

function expectedbackwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, riskmeasure::RiskMeasure, regularisation::Regularisation, backward_pass::BackwardPass)
    for t=(T-1):-1:1
        for pass = 1:getn(m)
            setrhs!(m, pass, t)
            for i=1:M
                # set_nonregularised_objective!(regularisation, X, subproblem(m, t+1, i))
                setexpectedvaluerhs!(m, subproblem(m, t+1, i))
                backsolve!(subproblem(m, t+1, i), 1)
            end
            xbar = getx(m, pass, t)
            markov = getmarkov(m, pass, t)
            cut = Cut(length(xbar))
            for i=1:M
                sd = stagedata(m, t+1, i)
                cut.intercept += transitionprobability(m, t, markov, i) * sd.objective_values[1]
                for v = 1:length(xbar)
                    cut.intercept -= transitionprobability(m, t, markov, i) * sd.dual_values[1][v] * xbar[v]
                    cut.coefficients[v] += transitionprobability(m, t, markov, i) * sd.dual_values[1][v]
                end
            end
            add_cut!(X, stagecut(m, t, markov), cut)
            addcut!(X, subproblem(m, t, markov), cut)
        end
    end
end

function expectedbackwardpass!{T, M, S, X, TM}(log, m::SDDPModel{T, M, S, X, TM}, riskmeasure::RiskMeasure, regularisation::Regularisation, backward_pass::BackwardPass)
    tic()
    expectedbackwardpass!(m, riskmeasure, regularisation, backward_pass)
    log.cuts += 1
    log.time_backwards += toq()
end
