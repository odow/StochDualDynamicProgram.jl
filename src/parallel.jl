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
            m.valuetogobound
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
function distribute_work!(results, f::Function, args...)
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                results[i] = remotecall_fetch(f, procid, args...)
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
function worker_simulate!(stagecuts, n::Int, vars::Vector{Symbol})
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
    nn = ceil(Int, n / length(workers()))
    distribute_work!(results, worker_simulate!, m.stagecuts, nn, vars)
    reduce_simulation!(m, results)
end

# ------------------------------------------------------------------------------
#
#  Monte Carlo estimation
#
function workermontecarloestimation(stagecuts, n, antitheticvariates)
    updatecuts!(stagecuts)
    montecarloestimation(antitheticvariates, m, n)
end

function parallelmontecarloestimation{T, M, S, X, TM}(antitheticvariates, m::SDDPModel{T, M, S, X, TM}, n::Int)
    results = Array(Vector{Float64}, length(workers()))
    nn = ceil(Int, n / length(workers()))
    distribute_work!(results, workermontecarloestimation, m.stagecuts, nn, antitheticvariates)
    return vcat(results...)
end


# ------------------------------------------------------------------------------
#
#  Forward Pass
#
function workerforwardpass!(stagecuts, n::Int, cutselection::CutSelectionMethod, forwardpass::ForwardPass)
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
    nn = ceil(Int, n / length(workers()))
    distribute_work!(results, workerforwardpass!, m.stagecuts, nn, cutselection, forwardpass)
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

function solvestage!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t, N, riskmeasure)
    # println("solving stage on processor $(myid())")
    tuples = [(pass, i, s) for i=1:M, s=1:S, pass=1:N]
    cuts = pmap(solvescenario!, repeated(t), tuples[:])
    cutsout = Array(Cut, N)
    x = 1:Int(M*S)
    for pass=1:N
        cutsout[pass] = reducecutvectors(m, cuts[(pass-1)*length(x) + x], riskmeasure, getmarkov(m, pass, t), t)
    end
    return cutsout
end

function reducecutvectors{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, x::Vector, riskmeasure::RiskMeasure, i, t)
    intercept = 0.
    coeffs = zeros(length(x[1].coefficients))
    idx=1
    for j=1:M
        for s=1:S
            stagedata(m, t, i).weightings_matrix[idx] = transitionprobability(m, t, i, j)*m.scenario_probability[s]
            idx+=1
        end
    end
    nestedcvar!(stagedata(m, t, i).weightings_matrix[:], Float64[c.intercept for c in x], riskmeasure.beta, riskmeasure.lambda, isa(X, Type{Max}))
    for cutidx in 1:length(x)
        intercept += x[cutidx].intercept * stagedata(m, t, i).weightings_matrix[cutidx]
        coeffs += x[cutidx].coefficients * stagedata(m, t, i).weightings_matrix[cutidx]
    end
    Cut(intercept / length(x), coeffs ./ length(x))
end

function solvescenario!(t, tuparg)
    pass, i, s = tuparg
    setrhs!(m, pass, t, i)
    load_scenario!(subproblem(m,t+1,i), s)
    backsolve!(subproblem(m,t+1,i), s)
    sd       = stagedata(m, t+1, i)
    thetahat = sd.objective_values[s]
    pihat    = [d[s] for d in sd.dual_values]
    xbar     = getx(m, pass, t)
    # println("(t, pass, i, s) = ($t, $pass, $i, $s): $thetahat")
    return Cut(thetahat - dot(xbar, pihat), pihat)
end
