# Copyright Oscar Dowson, 2017

function solvesubproblem!(pi::Vector{Float64}, m::JuMP.Model)
    status = solve(m)
    if status != :Optimal
        error("Subproblems must be solved to optimality")
    end
    for i = 1:length(ext(m).states)
        pi[i] = getdual(ext(m).states)
    end
    return getobjectivevalue(m)
end

function solvesubproblem(m::JuMP.Model)
    pi = zeros(length(ext(m).states))
    obj = solvesubproblem!(pi, m)
    return obj, pi
end

function setrhs!(m::JuMP.Model, rhs::Scenario)
    @assert length(rhs.con == length(rhs.values)
    @inbounds for i=1:length(rhs.con)
        JuMP.setrhs!(rhs.con[i], rhs.values[i])
    end
end

function convexcombination(r1::Rib, r2::Rib, price::Float64)
    lambda = (r2.price - price) / (r2.price - r1.price)
    return lambda * r1.theta + (1-lambda) * r2.theta
end

function interpolateribs(x::Vector{Rib}, price::Float64)
    if price < x[1].price
        warn("Price $(price) is outside discretisation. Extrapolating.")
        return convexcombination(x[1], x[2], price)
    end
    for i in 2:length(x)
        if price <= x[i].price
            return convexcombination(x[i-1], x[i], price)
        end
    end
    warn("Price $(price) is outside discretisation. Extrapolating.")
    return convexcombination(x[end-1], x[end], price)
end

function setprice!(m::JuMP.Model, scenario::PriceScenario, old_price::Float64)
    # step price forward
    newprice = scenario.dynamics(old_price, scenario.noise)
    # stage profit + convex combination of value to go functions
    _setobjective!(m, scenario.objective(newprice) + interpolateribs(ext(m).theta, newprice))
end

function addcut!{S,C,R}(m::SDDPModel{S,C,R}, t::Int, i::Int, p::Int, storage_idx)
    base_sp = stageproblem(m, t, i)
    # get the cut from the cut generator
    cut = cutgenerator(m.riskmeasure, S, getstate(base_sp), view(m.storage.pi, 1:storage_idx), view(m.storage.obj, 1:storage_idx), view(m.storage.probability, 1:storage_idx), t, i)

    # save the cut in the oracle
    storecut!(m.cutoracle, m, t, i, p, cut)

    # add the constraint
    addcutconstraint!(base_sp, cut, p)
end

function addcutconstraint!(m::JuMP.Model, cut::Cut, p::Int)
    ex = ext(m)
    if getobjectivesense(m) == :Min
        @constraint(m, ex.ribs[p].theta >= cut.intercept + dot(ex.states, cut.coefficients))
    else
        @constraint(m, ex.ribs[p].theta <= cut.intercept + dot(ex.states, cut.coefficients))
    end
end

function solvestage{S,C,R}(m::SDDPModel{S,C,R}, t::Int, i::Int, p::Int, x)
    base_sp = stageproblem(m, t, i)
    base_ex = ext(base_sp)
    markov_probability   = 0.0
    price_probability    = 0.0
    scenario_probability = 0.0
    storage_idx = 0
    # for each markov state
    for j in 1:nummarkovstates(m, t+1)
        # begin building probability
        markov_probability = transition(m.transition, t, i, j)
        sp = stageproblem(m, t+1, j)
        ex = ext(sp)
        for pnew in 1:length(ex.pricescenarios)
            price_probability = ex.pricescenarios[pnew].probability

            setprice!(sp, ex.pricescenarios[pnew], base_ex.ribs[p])
            for s in 1:length(ex.scenarios)
                scenario_probability = ex.scenarios[s].probability
                setrhs!(sp, ex.scenarios[s])

                storage_idx += 1
                obj = solvesubproblem!(m.storage.pi[storage_idx], sp)
                m.storage.obj[storage_idx] = obj
                m.storage.probability[storage_idx] = markov_probability * price_probability * scenario_probability
            end
        end
    end

    addcut!(m, t, i, p, storage_idx)

    # This can only be done if the price scenarios do not depend on markov state
    # delta::Float64
    # for ii in 1:nummarkovstates(m, t)
    #     i == ii && continue # skip current state
    #     k = 0
    #     for j in 1:nummarkovstates(m, t+1)
    #         delta = transition(m, t, ii, j) / transition(m, t, i, j)
    #         for pnew in 1:length(ex.pricescenarios)
    #             for s in 1:length(ex.scenarios)
    #                 k += 1
    #                 m.storage.probability[k] *= delta
    #             end
    #         end
    #     end
    #     addcut!(m, t, ii, p, storage_idx)
    # end
end

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
function backwardpass!(m::SDDPModel, markov_states::Vector{Int}, price_states::Vector{Int})
#
#     # solve final stage markov problems
#     for sp in m.stageproblems[end]
#         # solve markov state
#     end
#
end
