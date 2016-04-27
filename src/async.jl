# problem: As the number of cuts gets bigger, the number of processors involves decreases
# since the master process takes longer to add the cuts than it takes to recompute N more
# This means that eventually, it farms out to 1 process at a time while it adds those cuts

function async_worker_backward_pass!{T}(sc::Array{T,2}, n::Int)
    println("Proc $(myid()) received $(length(sc[1,1].cuts)) cuts")
    m.stagecuts = deepcopy(sc)
    oldn = Array(Int, (m.stages, m.markov_states,2))
    for stage=1:m.stages
        for markovstate=1:m.markov_states
            oldn[stage, markovstate, 1] = length(m.stagecuts[stage,markovstate].cuts)
            oldn[stage, markovstate, 2] = length(m.stagecuts[stage,markovstate].samplepoints)
        end
    end
    rebuild_stageproblems!(m)
    obj = zeros(n)
    for i=1:n
        obj[i] = backward_pass!(m, true)
    end
    N = getN(m.stagecuts[1,1])
    res = Array(Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, (m.stages, m.markov_states))
    for stage=1:m.stages
        for markovstate=1:m.markov_states
            res[stage, markovstate] = gettuple(m.stagecuts[stage, markovstate], oldn[stage,markovstate,1], oldn[stage,markovstate,2])
        end
    end
    res, obj
end
getN{N}(x::StageCuts{N}) = N
gettuple{N}(sc::StageCuts{N}, oldcuts::Int, oldsamplepoints::Int) = (sc.samplepoints[(oldsamplepoints+1):end], sc.cuts[(oldcuts+1):end])

function merge_stagecuts!{N}(m::SDDPModel, s::StageCuts{N}, stagecuts::Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}})
    s.samplepoints = unique(union(s.samplepoints, stagecuts[1]))
    s.cuts = unique(union(s.cuts, stagecuts[2]))
    recalculate_dominance!(m.sense, s)
end

function reduce_async_backwards_pass!{N}(m::SDDPModel, results::Array{Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, 2})
    for stage=1:m.stages
        for markovstate=1:m.markov_states
            merge_stagecuts!(m, m.stagecuts[stage, markovstate], results[stage, markovstate])
        end
    end
    m
end

function update_model!(m::SDDPModel, results, obj)
    reduce_async_backwards_pass!(m, results)
    rebuild_stageproblems!(m)

    # Recalculate bound
    solve_all_stage_problems!(m, 1)
    set_valid_bound!(m)

    test_and_set_ci!(m, obj)
end

function async_solve{N1<:Real, N2<:Real}(mm::SDDPModel, simulation_passes::Int, convergence_test_frequency::Int, maximum_iterations::Int, beta_quantile::N1, risk_lambda::N2,  cut_selection_frequency::Int, convergence_test::Bool, cuts_per_processor::Int)
    # Intialise model on worker processors
    nworkers = initialise_workers!(mm)

    # Initialise cut selection storage since we need it to communicate between processors
    initialise_cut_selection!(mm)

    # Set risk aversion parameters
    mm.beta_quantile, mm.risk_lambda = beta_quantile, risk_lambda

    total_iterations = 0
    objectives = RollingStorage(convergence_test_frequency)

    function updatemodel!(results, obj)
        total_iterations+=cuts_per_processor

        # Update simulated bound
        insert!(objectives, obj)

        # Combine cuts into model
        update_model!(mm, results, objectives)

        # Output to user
        # print_stats(m, objectives.n)
        println("$(total_iterations), $(mm.valid_bound), $(mm.confidence_interval)")

        return
    end
    get_stagecuts() = (mm.stagecuts)
    function check_continue()
        total_iterations < maximum_iterations && !terminate(rtol(mm) < 0., convergence_test)
    end
    @sync begin
        for procid in workers()
            @async begin
                while check_continue()
                    stagecuts = get_stagecuts()
                    results, obj = remotecall_fetch(async_worker_backward_pass!, procid, stagecuts, cuts_per_processor)
                    updatemodel!(results, obj)
                end
            end
        end
    end

    if terminate(rtol(mm) < 0., convergence_test)
        return CONVERGENCE_TERMINATION
    else
        return ITERATION_TERMINATION
    end
end

type RollingStorage{T,N}
    values::Vector{T}   # values
    i::Int              # next index
    n::Int              # num values
    sum::Float64        # mean of stored values
    sumsquares::Float64 # sum of squares of stored values
end
RollingStorage(T::DataType, n::Int) = RollingStorage{T,n}(zeros(T, n), 1, 0, 0., 0.)
RollingStorage(n::Int) = RollingStorage(Float64, n)

function Base.insert!{T,N}(s::RollingStorage{T,N}, x::T)
    s.sum += (x - s.values[s.i])            # Update sum
    s.sumsquares += (x^2- s.values[s.i]^2)  # Update sum of squares
    s.values[s.i] = x                       # Store value
    s.i += 1                                # Increment index
    if s.i > N
        s.i = 1                      # Wrap around if necessary
    end
    if s.n < N
        s.n += 1                     # Update total values stored
    end
    return
end
function Base.insert!{T,N}(s::RollingStorage{T,N}, X::Vector{T})
    for x in X
        insert!(s, x)
    end
end
Base.mean{T,N}(s::RollingStorage{T,N}) = s.sum / s.n
function Base.std{T,N}(s::RollingStorage{T,N})
    μ = mean(s)
    sqrt((s.sumsquares - 2 * μ * s.sum + s.n * μ^2) / (s.n-1))
end
Base.length{T,N}(s::RollingStorage{T,N}) = s.n
beta_quantile{T,N}(s::RollingStorage{T,N}, beta) = beta_quantile(s.values[1:s.n], beta)
