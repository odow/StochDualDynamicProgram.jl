function parallel_simulate!(m::SDDPModel, n::Int, vars::Vector{Symbol}=Symbol[])
    results = Array(Dict{Symbol, Any}, length(workers()))
    nn = ceil(Int, n / length(workers()))
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                results[i] = remotecall_fetch(worker_simulate!, procid, m.stagecuts, nn, vars)
            end
        end
    end
    reduce_simulation!(m, results)
end

function worker_simulate!{N}(sc::StageCuts{N}, n::Int, vars::Vector{Symbol}=Symbol[])
    m.stagecuts = deepcopy(sc)
    rebuild_stageproblems!(m)
    simulate(m, n, vars)
end

function reduce_simulation!(m::SDDPModel, results::Vector{Dict{Symbol, Any}})
    for res in results[2:end]
        merge_dicts!(m, results[1], res)
    end
    return results[1]
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


function parallel_forward_pass!(m::SDDPModel, npasses::Int=1)
    results = Array(Array{Float64}, length(workers()))
    nn = ceil(Int, npasses / length(workers()))
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                results[i] = remotecall_fetch(worker_forward_pass!, procid, m.stagecuts, nn)
            end
        end
    end
    # set new lower bound
    test_and_set_ci!(m, vcat(results...))
    return (true, npasses)
end

function worker_forward_pass!{T}(sc::Array{T,2}, n::Int)
    m.stagecuts = deepcopy(sc)
    forward_pass_kernel!(m, n)
end
