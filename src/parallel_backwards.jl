# function parallel_backward_pass!(m::SDDPModel, n::Int)
#     results = Array(Array{StageCuts{N}, 2}, length(workers()))
#     @sync begin
#         for (i, procid) in enumerate(workers())
#             @async begin
#                 results[i] = remotecall_fetch(worker_backward_pass!, procid, n)
#             end
#         end
#     end
#     reduce!(m, results)
# end
getN{N}(x::StageCuts{N}) = N
function parallel_backward_pass!(m::SDDPModel, n::Int)
    N = getN(m.stagecuts[1,1])
    results = Array(Array{Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, 2}, length(workers()))
    nn = ceil(Int, n / length(workers()))
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                results[i] = remotecall_fetch(worker_backward_pass!, procid, m.stagecuts, nn)
            end
        end
    end
    reduce!(m, results)
    rebuild_stageproblems!(m)
    solve_all_stage_problems!(m, 1)
    set_valid_bound!(m)
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
end

function sendtoworkers!(;kwargs...)
    for procid in workers()
        for (key, val) in kwargs
            @spawnat(procid, eval(StochDualDynamicProgram, Expr(:(=), key, val)))
        end
    end
end

# function worker_backward_pass!{N}(sc::Array{StageCuts{N},2}, n::Int)
#     m.stagecuts = sc
#     rebuild_stageproblems!(m)
#     for i=1:n
#         backward_pass!(m, true)
#     end
#     m.stagecuts
# end
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
# gettuple{N}(sc::StageCuts{N}, oldcuts::Int, oldsamplepoints::Int) = (sc.samplepoints[min(length(sc.samplepoints),oldsamplepoints+1):end], sc.cuts[min(length(sc.cuts,oldcuts+1):end])

# function reduce!{N}(m::SDDPModel, results::Vector{Array{StageCuts{N}, 2}})
#     for stage=1:m.stages
#         for markovstate=1:m.markov_states
#             merge!(m, m.stagecuts[stage, markovstate], [res[stage, markovstate] for res in results])
#         end
#     end
#
# end
function reduce!{N}(m::SDDPModel, results::Vector{Array{Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, 2}})
    for stage=1:m.stages
        for markovstate=1:m.markov_states
            merge!(m, m.stagecuts[stage, markovstate], [res[stage, markovstate] for res in results])
        end
    end

end

# function merge!{N}(m::SDDPModel, s::StageCuts{N}, stagecuts::Vector{StageCuts{N}})
#     s.samplepoints = unique(union([sc.samplepoints for sc in stagecuts]...))
#     s.cuts = unique(union([sc.cuts for sc in stagecuts]...))
#     s.n = length(s.samplepoints)
#     s.nondominated = zeros(length(s.cuts))
#     s.activecut = zeros(s.n)
#
#     for (sample_idx, samplepoint) in enumerate(s.samplepoints)
#         best_idx = 1
#         best_y = evaluate(s.cuts[1], samplepoint)
#         for cut_idx in 2:length(s.cuts)
#             y = evaluate(s.cut[cut_idx], samplepoint)
#             if is_dominated(m.sense, best_y, y)
#                 best_idx = cut_idx
#                 best_y = y
#             end
#         end
#         s.nondominated[best_idx] += 1
#         s.activecut[sample_idx] = best_idx
#     end
# end
function merge!{N}(m::SDDPModel, s::StageCuts{N}, stagecuts::Vector{Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}})
    s.samplepoints = unique(union(s.samplepoints, [sc[1] for sc in stagecuts]...))
    s.cuts = unique(union(s.cuts, [sc[2] for sc in stagecuts]...))
    s.n = length(s.samplepoints)
    s.nondominated = zeros(length(s.cuts))
    s.activecut = zeros(s.n)

    for (sample_idx, samplepoint) in enumerate(s.samplepoints)
        best_idx = 1
        best_y = evaluate(s.cuts[1], samplepoint)
        for cut_idx in 2:length(s.cuts)
            y = evaluate(s.cuts[cut_idx], samplepoint)
            if is_dominated(m.sense, best_y, y)
                best_idx = cut_idx
                best_y = y
            end
        end
        s.nondominated[best_idx] += 1
        s.activecut[sample_idx] = best_idx
    end
end
