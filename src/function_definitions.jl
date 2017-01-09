# Copyright Oscar Dowson, 2017

dominates(::Maximisation, x, y) = x < y
dominates(::Minimisation, x, y) = x > y

# forward pass
#     with standard probabilities
#     with uniform probabilities
#
#     scenario incrementation
#
#
# backward pass
#     sample trajectory
#     solve all subproblems
#     add cut
#     back track
#
# cut selection
#
# simulation
#
# function backwardpass!(m::SDDPModel)
#
#     # solve final stage markov problems
#     for sp in m.stageproblems[end]
#         # solve markov state
#     end
#
# end
function setrhs!(m::JuMP.Model, rhs::Scenario)
    @assert rhs.con == rhs.values
    for i=1:length(rhs.con)
        JuMP.setrhs!(rhs.con[i], rhs.values[i])
    end
end
function setprice!(m::JuMP.Model, price, noise)
    newprice = ext(m).pricescenarios.dynamics(price, noise)
    ext(m).pricescenarios.objective(newprice)
end

function solvestage(m::SDDPModel, t::Int, i::Int, p::Int, x)
    sp = stageproblem(m, t, i)
    ex = ext(sp)
    for rhs in ex.scenarios.samples
        setrhs!(sp, rhs)
        # for price in ex.pricescenarios

        # end
    end
    # for each rib
    #     for each markov state
    #         for each scenario
    #             solve
    #             cache theta, π
    #         end
    #     end
    #     if pricescenario independent of markov state? # can check at beginning
    #         calculate cuts across markov states
    #     else
    #         just add cut to this markov state
    #     end
    # end
end

function solvesubproblem!(pi::Vector{Float64}, m::JuMP.Model)
    status = solve(m)
    if status != :Optimal
        error("Subproblems must be solved to optimality")
    end
    map!(JuMP.getdual, pi, ext(m).states)
    return getobjectivevalue(m)
end
function solvesubproblem(m::JuMP.Model)
    pi = zeros(length(ext(m).states))
    obj = solvesubproblem!(pi, m)
    return obj, pi
end
