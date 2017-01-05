# Copyright Oscar Dowson, 2017

dominates(::Maximisation, x, y) = x < y
dominates(::Minimisation, x, y) = x > y

forward pass
    with standard probabilities
    with uniform probabilities

    scenario incrementation


backward pass
    sample trajectory
    solve all subproblems
    add cut
    back track

cut selection

simulation

function backwardpass!(m::SDDPModel)

    # solve final stage markov problems
    for sp in m.stageproblems[end]
        # solve markov state
    end

end

function solvestage(m::SDDPModel, t::Int, x::incomingstate)
    for each rib
        for each markov state
            for each scenario
                solve
                cache theta, π
            end
        end
        if pricescenario independent of markov state? # can check at beginning
            calculate cuts across markov states
        else
            just add cut to this markov state
        end
    end
end
