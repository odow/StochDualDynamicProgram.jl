function sendtoworkers!(;kwargs...)
    for procid in workers()
        for (key, val) in kwargs
            @spawnat(procid, eval(StochDualDynamicProgram, Expr(:(=), key, val)))
        end
    end
end

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
    return length(workers())
end

function distribute_work!(results, f::Function, args...)
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                results[i] = remotecall_fetch(f, procid, args...)
            end
        end
    end
end
function distribute_work_void!(f::Function, args...)
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                remotecall_fetch(f, procid, args...)
            end
        end
    end
end

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
    y = Array(Vector{Float64}, (T,M))
    for t=1:T
        for i=1:M
            y[t,i] = stagecut(m, t, i).samplepoints[(end-n+1):end]
        end
    end
    return y
end
getlatestsamplepoints{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, cutselection::CutSelectionMethod, n) = Array(Vector{Float64}, (0,0))

function reduceforwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, results::Vector{Tuple{ForwardPassData, Array{Vector{Float64}, 2}}})
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
                    push!(stagecut(m, t, i).samplepoints, result[2][t,i]...)
                end
            end
        end
    end
end

function parallelforwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int, cutselection::CutSelectionMethod, forwardpass::ForwardPass)
    results = Array(Tuple{ForwardPassData, Array{Vector{Float64}, 2}}, length(workers()))
    nn = ceil(Int, n / length(workers()))
    distribute_work!(results, workerforwardpass!, m.stagecuts, nn, cutselection, forwardpass)
    reduceforwardpass!(m, results)
end


# ------------------------------------------------------------------------------
#
#  Backward Pass
#
function workerbackwardpass!(pass, T, riskmeasure, regularisation)
    cuts = Array(Cut, T-1)
    for t=(T-1):-1:1
        setrhs!(m, pass, t)
        solveall!(m, t+1, regularisation)
        reweightscenarios!(m, t, getmarkov(m, pass, t), riskmeasure.beta, riskmeasure.lambda)
        cuts[t] = addcut!(m, pass, t, getmarkov(m, pass, t))
    end
    return cuts
end

function reducebackwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, results)
    for pass=1:getn(m)
        for t=1:(T-1)
            # add to cutselection storage
            add_cut!(X, stagecut(m, t, getmarkov(m, pass, t)), results[pass][t])
            # add to problem
            addcut!(X, subproblem(m, t, getmarkov(m, pass, t)), results[pass][t])
        end
    end
end

function parallelbackwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, riskmeasure::RiskMeasure, regularisation::Regularisation=NoRegularisation())
    updateforwardstorage!(m.forwardstorage)
    # This returns a vector of vector of cuts
    results = pmap(workerbackwardpass!, 1:getn(m), repeated(T), repeated(riskmeasure), repeated(regularisation))
    reducebackwardpass!(m, results)
end
