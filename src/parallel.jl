function sendtoworkers!(;kwargs...)
    for procid in workers()
        for (key, val) in kwargs
            @spawnat(procid, eval(StochDualDynamicProgram, Expr(:(=), key, val)))
        end
    end
end

function initialise_workers!(m::SDDPModel)
    sendtoworkers!(m = SDDPModel(
        sense=isa(m.sense, Val{:Min})?(:Min):(:Max),
        stages=m.stages,
        markov_states=m.markov_states,
        scenarios=m.scenarios,
        scenario_probability=deepcopy(m.scenario_probability.values),
        transition=deepcopy(m.transition),
        initial_markov_state=m.init_markov_state,
        conf_level=m.QUANTILE,
        solver=m.LPSOLVER,
        value_to_go_bound=m.value_to_go_bound,
        cuts_filename=nothing,
        build_function=m.build_function!
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

# ------------------------------------------------------------------------------
#
#   Backwards pass functionality
#
function worker_backward_pass!{T}(sc::Array{T,2}, n::Int)
    m.stagecuts = deepcopy(sc)
    oldn = Array(Int, (m.stages, m.markov_states,2))
    for stage=1:m.stages
        for markovstate=1:m.markov_states
            oldn[stage, markovstate, 1] = length(m.stagecuts[stage,markovstate].cuts)
            oldn[stage, markovstate, 2] = length(m.stagecuts[stage,markovstate].samplepoints)
        end
    end
    rebuild_stageproblems!(m)
    for i=1:n
        backward_pass!(m, true)
    end
    N = getN(m.stagecuts[1,1])
    res = Array(Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, (m.stages, m.markov_states))
    for stage=1:m.stages
        for markovstate=1:m.markov_states
            res[stage, markovstate] = gettuple(m.stagecuts[stage, markovstate], oldn[stage,markovstate,1], oldn[stage,markovstate,2])
        end
    end
    res
end
gettuple{N}(sc::StageCuts{N}, oldcuts::Int, oldsamplepoints::Int) = (sc.samplepoints[(oldsamplepoints+1):end], sc.cuts[(oldcuts+1):end])

function merge_stagecuts!{N}(m::SDDPModel, s::StageCuts{N}, stagecuts::Vector{Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}})
    s.samplepoints = unique(union(s.samplepoints, [sc[1] for sc in stagecuts]...))
    s.cuts = unique(union(s.cuts, [sc[2] for sc in stagecuts]...))
    recalculate_dominance!(m.sense, s)
end

function reduce_backwards_pass!{N}(m::SDDPModel, results::Vector{Array{Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, 2}})
    for stage=1:m.stages
        for markovstate=1:m.markov_states
            merge_stagecuts!(m, m.stagecuts[stage, markovstate], [res[stage, markovstate] for res in results])
        end
    end
end

getN{N}(x::StageCuts{N}) = N
function parallel_backward_pass!(m::SDDPModel, n::Int)
    N = getN(m.stagecuts[1,1])
    results = Array(Array{Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, 2}, length(workers()))
    nn = ceil(Int, n / length(workers()))
    distribute_work!(results, worker_backward_pass!, m.stagecuts, nn)
    reduce_backwards_pass!(m, results)
    rebuild_stageproblems!(m)
    solve_all_stage_problems!(m, 1)
    set_valid_bound!(m)
end

# ------------------------------------------------------------------------------
#
#   Forwards pass functionality
#
function worker_forward_pass!{T}(sc::Array{T,2}, n::Int)
    m.stagecuts = deepcopy(sc)
    forward_pass_kernel!(m, n)
end

function parallel_forward_pass!(m::SDDPModel, npasses::Int=1)
    results = Array(Array{Float64}, length(workers()))
    nn = ceil(Int, npasses / length(workers()))
    distribute_work!(results, worker_forward_pass!, m.stagecuts, nn)
    # set new lower bound
    test_and_set_ci!(m, vcat(results...))
    return (rtol(m) < 0., npasses)
end

# ------------------------------------------------------------------------------
#
#   Simulation functionality
#
function worker_simulate!{T}(sc::Array{T,2}, n::Int, vars::Vector{Symbol}=Symbol[])
    m.stagecuts = deepcopy(sc)
    rebuild_stageproblems!(m)
    simulate(m, n, vars)
end

function merge_dicts!(m::SDDPModel, d1::Dict{Symbol, Any}, d2::Dict{Symbol, Any})
    for key in keys(d1)
        @assert haskey(d2, key)
        if key == :Objective
            d1[key] = vcat(d1[key], d2[key])
        else
            for stage=1:m.stages
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

function parallel_simulate(m::SDDPModel, n::Int, vars::Vector{Symbol}=Symbol[])
    results = Array(Dict{Symbol, Any}, length(workers()))
    nn = ceil(Int, n / length(workers()))
    distribute_work!(results, worker_simulate!, m.stagecuts, nn, vars)
    reduce_simulation!(m, results)
end
