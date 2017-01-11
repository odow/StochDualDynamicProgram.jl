# Copyright Oscar Dowson, 2017

function solvetooptimality!(m::JuMP.Model)
    status = solve(m)
    if status != :Optimal
        error("Subproblems must be solved to optimality")
    end
end

function solvesubproblem!(pi::Vector{Float64}, m::JuMP.Model)
    solvetooptimality!(m)
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
    @assert length(rhs.con) == length(rhs.values)
    @inbounds for i=1:length(rhs.con)
        JuMP.setrhs!(rhs.con[i], rhs.values[i])
    end
end

function setstate!(m::JuMP.Model, state::Vector{Float64})
    ex = ext(m)
    for i = 1:length(state)
        JuMP.setrhs!(ex.states[i].con, state[i])
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

    newprice
end
setprice!(m::JuMP.Model, scenario::PriceScenario, old_rib::Rib) = setprice!(m, scenario, old_rib.price)

function addcut!{S,C,R}(m::SDDPModel{S,C,R}, t::Int, i::Int, p::Int, storage_idx)
    base_sp = stageproblem(m, t, i)
    # get the cut from the cut generator
    cut = cutgenerator(m.riskmeasure, S, getstate(base_sp), view(m.storage.pi, 1:storage_idx), view(m.storage.obj, 1:storage_idx), view(m.storage.probability, 1:storage_idx), t, i)

    # save the cut in the oracle
    storecut!(m.cutoracle, m, t, i, p, cut)

    # add the constraint
    addcutconstraint!(base_sp, cut, p)

    return cut
end

function addcutconstraint!(m::JuMP.Model, cut::Cut, p::Int)
    ex = ext(m)
    if getobjectivesense(m) == :Min
        @constraint(m, ex.ribs[p].theta >= cut.intercept + dot(ex.states, cut.coefficients))
    else
        @constraint(m, ex.ribs[p].theta <= cut.intercept + dot(ex.states, cut.coefficients))
    end
end

function solvestage!{S,C,R}(m::SDDPModel{S,C,R}, cutstore::Vector{CutContainer}, t::Int, i::Int, p::Int, x)
    base_sp = stageproblem(m, t, i)
    markov_probability   = 0.0
    price_probability    = 0.0
    scenario_probability = 0.0
    storage_idx = 0
    # for each markov state
    for j in 1:nummarkovstates(m, t+1)
        # begin building probability
        markov_probability = transition(m, t, i, j)
        sp = stageproblem(m, t+1, j)
        setstate!(sp, x)

        for pnew in 1:numpricescenarios(sp)
            price_probability = probability(pricescenario(sp, pnew))
            setprice!(sp, pricescenario(sp, pnew), rib(base_sp, p))
            for s in 1:numscenarios(sp)
                scenario_probability = probability(scenario(sp, s))
                setrhs!(sp, scenario(sp, s))
                storage_idx += 1
                obj = solvesubproblem!(dualstore(m, storage_idx), sp)
                storeobj!(m, storage_idx, obj)
                storeprob!(m, storage_idx, markov_probability * price_probability * scenario_probability)
            end
        end
    end
    cut = addcut!(m, t, i, p, storage_idx)
    push!(cutstore, CutContainer(cut, t, i, p))
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


function backwardpass!(m::SDDPModel, cutstore::Vector{CutContainer}, num_forward_passes::Int)
    # backtrack up stages
    for t in reverse(1:(numstages(m)-1))
        # for scenario incrementation
        for i=1:num_forward_passes
            markov = getmarkov(m, i, t)
            # it's critical to solve upper and lower rib
            for p in boundingribs(m, t, i)
                solvestage!(m, cutstore, t, markov, p, getstate(m, i, t))
            end
        end
    end
end
function boundingribs(m, t, pass)
    ex = ext(stageproblem(m, t, getmarkov(m, t, pass)))
    p = getprice(m, t, pass)
    if p < ex.theta[1].price
        return (1,2)
    elseif p > ex.theta[end].price
        return (length(ex.theta)-1, length(ex.theta))
    else
        for i=2:length(ex.theta)
            if p < ex.theta[i].price
                return (i-1, i)
            end
        end
    end
end


function samplefirstmarkov(m::SDDPModel)
    r = rand()
    for i=1:length(m.initialmarkovprobability)
        r -= m.initialmarkovprobability[i]
        if r <= 0.0
            return i
        end
    end
end
function getfirstprice(m::SDDPModel)
    if numpriceribs(m, 1) != 1
        @assert !isnan(m.firstprice)
    end
    return m.firstprice
end

function forwardpass!(m::SDDPModel, num_forward_passes::Int)
    for i=1:num_forward_passes
        if length(m.forwardstorage) < i
            push!(m.forwardstorage, ForwardStorage(m))
        end
        setmarkov!(m, i, 1, samplefirstmarkov(m))
        setprice!(m, i, 1, getfirstprice(m))
        solvestageforward!(m, 1, i)
        for t in 2:numstages(m)
            solvestageforward!(m, t, i)
        end
    end
end


function solvestageforward!(m::SDDPModel, stage, pass)
    # markov = samplemarkov(m.sampler, gettransitionmatrix(m, stage), stage, getmarkov(m, i, stage-1))
    # setmarkov!(m, i, t, markov) # sample markov

    sp = stageproblem(m, stage, getmarkov(m, pass, stage))
    setstate!(sp, getstate(m, pass, stage-1))

    setrhs!(sp, samplescenario(m.sampler, ext(sp).scenarios, stage, markov)) # sample rhs
    if numpricescenarios(sp) > 0
        newprice = setprice!(sp, samplepricescenario(m.sampler, ext(sp).pricescenarios, stage, markov), getprice(m, pass, stage-1)) # sample price
        setprice!(m, pass, stage, newprice)
    end
    solvetooptimality!(sp) # solve subproblem
    getstate!(sp, getstate(m, pass, stage))
    addstatevisited!(m, t, getstate(m, pass, stage))
end

function addcuts!(m::SDDPModel, cuts::Vector{CutContainer})
    for c in cuts
        addcut!(m, c)
    end
end
function addcut!(m::SDDPModel, cut::CutContainer)
        # save the cut in the oracle
        storecut!(m.cutoracle, m, c.stage, c.markov, c.rib, c.cut)
        # add the constraint
        addcutconstraint!(stageproblem(m, c.stage, c.markov), c.cut, c.rib)
    end
end


function asyncronousiteration!(m::SDDPModel, scenarios::Int, newcuts::Vector{CutContainer})
    addcuts!(m, newcuts) # add new cuts
    # cut selection
    rebuildsubproblems!(m)

    forwardpass!(m, scenarios)
    cutstore = CutContainer[]
    backwardpass!(m, cutstore, scenarios)
    # check termination
    terminate = false
    return terminate, cutstore
end
