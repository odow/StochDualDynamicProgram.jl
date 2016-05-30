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

function distribute_work!(results, f, args...)
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                results[i] = remotecall_fetch(f, procid, args...)
            end
        end
    end
end

function updatecuts!(stagecuts)
    m.stagecuts = deepcopy(stagecuts)
    rebuild_stageproblems!(m)
end
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
