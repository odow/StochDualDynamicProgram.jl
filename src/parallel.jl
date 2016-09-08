"""
    sendtoworkers!(;kwargs...)

This function spawns keyword variables on all the workers.

Usage:

    sendtoworkers!(a=1)
"""
function sendtoworkers!(;kwargs...)
    for procid in workers()
        for (key, val) in kwargs
            @spawnat(procid, eval(StochDualDynamicProgram, Expr(:(=), key, val)))
        end
    end
end

"""
    initialise_workers!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM})

Copies the SDDPModel `m` to all the workers.
"""
function initialise_workers!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM})
    sendtoworkers!(
        m = SDDPModel(
            m.build_function!,
            m.initial_markov_state,
            M,
            S,
            deepcopy(m.scenario_probability.values),
            isa(X, Min)?(:Min):(:Max),
            T,
            deepcopy(m.transition),
            m.solver,
            m.backupsolver,
            m.valuetogobound,
            m.roundingaccuracy
        )
    )
    distribute_work_void!(initialisecutstorage!)
    return length(workers())
end
initialisecutstorage!() = initialise_cutstorage!(m)

"""
    distribute_work!(results, f::Function, args...)

Runs the function `f` with the arguments `args` on all workers and stores the solution in `results`.
"""
function distribute_work!(results, f::Function, vecarg::Vector, args...)
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                results[i] = remotecall_fetch(f, procid, vecarg[i], args...)
            end
        end
    end
end

"""
    distribute_work_void!(f::Function, args...)

Runs the function `f` with the arguments `args` on all workers.
"""
function distribute_work_void!(f::Function, args...)
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                remotecall_fetch(f, procid, args...)
            end
        end
    end
end

"""
    updatecuts!(stagecuts)

Copies the stage cuts `stagecuts` to the model `m` on all workers (assumes `initialise_workers` has been run).
"""
function updatecuts!(stagecuts)
    m.stagecuts = deepcopy(stagecuts)
    rebuild_stageproblems!(m)
end

updateforwardstorageworker!(forwardstorage) = (m.forwardstorage = deepcopy(forwardstorage))
updateforwardstorage!(forwardstorage) = distribute_work_void!(updateforwardstorageworker!, forwardstorage)

# ------------------------------------------------------------------------------
#
#   Simulation functionality
#
function worker_simulate!(n::Int, stagecuts, vars::Vector{Symbol})
    updatecuts!(stagecuts)
    simulate(m, n, vars)
end

function merge_dicts!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, d1::Dict{Symbol, Any}, d2::Dict{Symbol, Any})
    for key in keys(d1)
        @assert haskey(d2, key)
        if key == :Objective
            d1[key] = vcat(d1[key], d2[key])
        else
            for stage=1:T
                d1[key][stage] = vcat(d1[key][stage], d2[key][stage])
            end
        end
    end
end

function reduce_simulation!(m::SDDPModel, results::Vector{Dict{Symbol, Any}})
    for res in results[2:end]
        merge_dicts!(m, results[1], res)
    end
    return results[1]
end

function parallelsimulate(m::SDDPModel, n::Int, vars::Vector{Symbol})
    nworkers = initialise_workers!(m)
    results = Array(Dict{Symbol, Any}, length(workers()))

    distribute_work!(results, worker_simulate!, divideup(n, length(workers())), m.stagecuts,vars)
    reduce_simulation!(m, results)
end

function divideup(n::Int, p::Int)
    xl = floor(Int, n / p)
    y = ones(Int, p) * xl
    for i=1:(n - sum(y))
        y[i] += 1
    end
    return y
end

# ------------------------------------------------------------------------------
#
#  Monte Carlo estimation
#
function workermontecarloestimation(n, stagecuts, antitheticvariates)
    updatecuts!(stagecuts)
    montecarloestimation(antitheticvariates, m, n)
end

function parallelmontecarloestimation{T, M, S, X, TM}(antitheticvariates, m::SDDPModel{T, M, S, X, TM}, n::Int)
    results = Array(Vector{Float64}, length(workers()))
    nn = ceil(Int, n / length(workers()))
    distribute_work!(results, workermontecarloestimation, divideup(n, length(workers())), m.stagecuts, antitheticvariates)
    return vcat(results...)
end


# ------------------------------------------------------------------------------
#
#  Forward Pass
#
function workerforwardpass!(n::Int, stagecuts, cutselection::CutSelectionMethod, forwardpass::ForwardPass)
    updatecuts!(stagecuts)
    force_resizeforwardstorage!(m, n) # resize storage for forward pass
    forwardpass!(m, n, cutselection, forwardpass)
    deepcopy(m.forwardstorage), getlatestsamplepoints(m, cutselection, n)
end

function getlatestsamplepoints{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, cutselection::LevelOne, n)
    y = Array(Vector{NTuple}, (T,M))
    for t=1:T
        for i=1:M
            y[t,i] = stagecut(m, t, i).samplepoints[max(1, end-n+1):end]
        end
    end
    return y
end
getlatestsamplepoints{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, cutselection::CutSelectionMethod, n) = Array(Vector{NTuple}, (0,0))

function reduceforwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, results::Vector{Tuple{ForwardPassData, Array{Vector{NTuple}, 2}}})
    bign = sum(map(r->r[1].n, results))
    force_resizeforwardstorage!(m, bign)
    pass = 1
    for result in results
        for littlen = 1:result[1].n
            m.forwardstorage.obj[pass] = result[1].obj[littlen]
            m.forwardstorage.W[pass]   = result[1].W[littlen]
            m.forwardstorage.x[pass]   = result[1].x[littlen]
            pass += 1
        end
        if size(result[2])[1] > 0
            for t=1:T
                for i=1:M
                    for res in result[2][t,i]
                        push!(stagecut(m, t, i).samplepoints, res)
                    end
                end
            end
        end
    end
end

function parallelforwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int, cutselection::CutSelectionMethod, forwardpass::ForwardPass)
    results = Array(Tuple{ForwardPassData, Array{Vector{NTuple}, 2}}, length(workers()))
    distribute_work!(results, workerforwardpass!, divideup(n, length(workers())), m.stagecuts, cutselection, forwardpass)
    reduceforwardpass!(m, results)
end


# ------------------------------------------------------------------------------
#
#  Backward Pass
#   solve trajectories
function reweightandgetcut!(m, t, i, pass, beta, lambda)
    reweightscenarios!(m, t,i, beta, lambda)
    addcut!(m, pass, t, i)
end

function setandsolve!(m, t, pass, regularisation)
    setrhs!(m, pass, t)
    solveall!(m, t+1, regularisation)
end

function addcuts!(X, m, t, i, result)
    add_cut!(X, stagecut(m, t, i), result)  # add to cutselection storage
    addcut!(X, subproblem(m, t, i), result) # add to problem
end

function workerbackwardpass!(pass, T, riskmeasure, regularisation)
    cuts = Array(Cut, T-1)
    for t=(T-1):-1:1
        setandsolve!(m, t, pass, regularisation)
        cuts[t] = reweightandgetcut!(m, t,  getmarkov(m, pass, t), pass, riskmeasure.beta, riskmeasure.lambda)
    end
    return cuts
end

function workerbackwardpassmulticut!(pass, T, M, riskmeasure, regularisation)
    cuts = Array(Cut, (T-1, M))
    for t=(T-1):-1:1
        setandsolve!(m, t, pass, regularisation)
        for i=1:M
            cuts[t, i] = reweightandgetcut!(m, t, i, pass, riskmeasure.beta, riskmeasure.lambda)
        end
    end
    return cuts
end

function reducebackwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, results)
    for t=1:(T-1)
        for pass=1:getn(m)
            addcuts!(X, m, t, getmarkov(m, pass, t), results[pass][t])
        end
    end
end

function reducebackwardpassmulticut!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, results)
    for t=1:(T-1)
        for pass=1:getn(m)
            for i=1:M
                addcuts!(X, m, t, i, results[pass][t, i])
            end
        end
    end
end

function parallelbackwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, riskmeasure::RiskMeasure, regularisation::Regularisation, backward_pass::BackwardPass)
    updateforwardstorage!(m.forwardstorage)
    # This returns a vector of vector of cuts
    if backward_pass.multicut
        results = pmap(workerbackwardpassmulticut!, 1:getn(m), repeated(T), repeated(M), repeated(riskmeasure), repeated(regularisation))
        reducebackwardpassmulticut!(m, results)
    else
        results = pmap(workerbackwardpass!, 1:getn(m), repeated(T), repeated(riskmeasure), repeated(regularisation))
        reducebackwardpass!(m, results)
    end
end

# solve all the t subproblems
function parallelsolveall!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t::Int, regularisation::Regularisation)
    pmap(workersolveall!,  [(i,s) for i=1:M, s=1:S][:], repeaded(X), repeated(t), repeated(regularisation))
end
function workersolveall!(i, X, t, regularisation)
    set_nonregularised_objective!(regularisation, X, subproblem(m,t,i[1]))
    load_scenario!(subproblem(m,t,i[1]), i[2])
    backsolve!(subproblem(m,t,i[1]), i[2])
end

# ------------------------------------------------------------------------------
#
#  Backward Pass Version 2.
#   solve subproblems
function bestparallelbackwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, riskmeasure::RiskMeasure, regularisation::Regularisation, backward_pass::BackwardPass)
    updateforwardstorage!(m.forwardstorage)
    set_nonregularised_objective_all!(m, regularisation, X, T, M)
    distribute_work_void!(set_nonregularised_objective_all!, regularisation, X, T, M)
    for t=(T-1):-1:1
        cuts = solvestage!(m, t, getn(m), riskmeasure)
        for pass=1:getn(m)
            addcuts!(X, m, t, getmarkov(m, pass, t), cuts[pass])
            distribute_work_void!(addremotecut!, X, t, getmarkov(m, pass, t), cuts[pass])
        end
    end
end
function set_nonregularised_objective_all!(m, regularisation, X, T, M)
    for t=1:T
        for i=1:M
            set_nonregularised_objective!(regularisation, X, subproblem(m, t, i))
        end
    end
end
set_nonregularised_objective_all!(regularisation, X, T, M) = set_nonregularised_objective_all!(m, regularisation, X, T, M)

addremotecut!(X, t, i, cut) = addcuts!(X, m, t, i, cut)

# function solvestage!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t, N, riskmeasure)
#     tuples = [(pass, i, s) for i=1:M, s=1:S, pass=1:N]
#     cuts = pmap(solvescenario!, repeated(t), tuples[:])
#     cutsout = Array(Cut, N)
#     x = 1:Int(M*S)
#     for pass=1:N
#         cutsout[pass] = reducecutvectors(m, cuts[(pass-1)*length(x) + x], riskmeasure, pass, t)
#     end
#     return cutsout
# end
function solvestage!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t, N, riskmeasure)
    tuples = [(pass, i) for i=1:M, pass=1:N]
    results = pmap(solvescenario!, repeated(t), tuples[:])
    cutsout = Array(Cut, N)
    for pass=1:N
        cutsout[pass] = reducecutvectors(m, results[((pass-1)*M+1):(pass*M)], riskmeasure, pass, t)
    end
    return cutsout
end

# function reducecutvectors{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, x::Vector, riskmeasure::RiskMeasure, pass, t)
#     cutout = Cut(0., zeros(length(x[1].coefficients)))
#     i = getmarkov(m, pass, t)
#     reweightscenarios!(m, Float64[c.intercept for c in x], t, i, riskmeasure.beta, riskmeasure.lambda)
#     for cutidx in 1:length(x)
#         prob = stagedata(m, t, i).weightings_matrix[cutidx]
#         cutout.intercept    += (x[cutidx].intercept - dot(x[cutidx].coefficients, getx(m, pass, t))) * prob
#         cutout.coefficients += x[cutidx].coefficients * prob
#     end
#     return cutout
# end
function reducecutvectors{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, results::Vector, riskmeasure::RiskMeasure, pass, t)
    objectives  = zeros(M*S)
    dual_values = Array(Vector{Float64}, M*S)
    idx = 0
    for s = 1:S
        for i = 1:M
            idx += 1
            objectives[idx]  = results[i][1][s]
            dual_values[idx] = results[i][2][s]
        end
    end
    cutout = Cut(0., zeros(length(dual_values[1])))
    i = getmarkov(m, pass, t)
    reweightscenarios!(m, objectives, t, i, riskmeasure.beta, riskmeasure.lambda)
    cutidx = 0
    for scenario=1:S
        for markov=1:M
            cutidx += 1
            prob = stagedata(m, t, i).weightings_matrix[markov, scenario]
            cutout.intercept    += (objectives[cutidx] - dot(dual_values[cutidx], getx(m, pass, t))) * prob
            cutout.coefficients += dual_values[cutidx] * prob
        end
    end
    cutout.coefficients = round(cutout.coefficients, m.roundingaccuracy)
    return cutout
end

# function solvescenario!(t, tuparg)
#     pass, i, s = tuparg
#     setrhs!(m, pass, t, i)
#     load_scenario!(subproblem(m,t+1,i), s)
#     backsolve!(subproblem(m,t+1,i), s)
#     sd       = stagedata(m, t+1, i)
#     # println("(t, pass, i, s) = ($t, $pass, $i, $s): $thetahat")
#     return Cut{length(sd.dual_values[1])}(sd.objective_values[s], sd.dual_values[s])
# end
function solvescenario!(t, tuparg)
    pass, i = tuparg
    setrhs!(m, pass, t, i)
    solvescenarios!(m, t+1, i)
    return copy(stagedata(m, t+1, i).objective_values), deepcopy(stagedata(m, t+1, i).dual_values)
end
