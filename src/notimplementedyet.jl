macro newstate(sp, x, x0)
    @assert x0.head == :kw # must be a keyword
    symin = x0.args[1]     # name of the statein variable
    rhs   = x0.args[2]     # values for the statein variable
    if Base.Meta.isexpr(x, :comparison) # if its a comparison
        if length(x.args) == 5          # double sided
            xin = x.args[3]             # variable is in middle
        elseif length(x.args) == 3      # single comparison
            xin = x.args[1]             # variable is on left
        else
            error("Unknown format for $(x)")
        end
    else
        xin = x                         # no bounds
    end
    quote
        # create unbounded statein variables
        if isa($xin, Expr)
            $(xin).args[1] = $symin
        end
        statein = @variable($(esc(sp)), $(esc(xin)))
        # create jump variables and initialise with starting solution
        stateout = @variable($(esc(sp)), $(esc(x)), start=$(esc(rhs)))
        # add to stage data
        if typeof(stateout) != JuMP.Variable
            push!(stagedata($(esc(sp))).state_vars, stateout[:]...)
            for i=1:length(stateout[:])
                push!(stagedata($(esc(sp))).dual_constraints, @constraint($(esc(sp)), statein[i] == getvalue(stateout[i])))
            end
        else # single state variable
            push!(stagedata($(esc(sp))).state_vars, stateout)
            push!(stagedata($(esc(sp))).dual_constraints, @constraint($(esc(sp)), statein == getvalue(stateout)))
        end
    end
end

# # problem: As the number of cuts gets bigger, the number of processors involves decreases
# # since the master process takes longer to add the cuts than it takes to recompute N more
# # This means that eventually, it farms out to 1 process at a time while it adds those cuts
#
# function async_worker_backward_pass!{T}(sc::Array{T,2}, n::Int)
#     println("Proc $(myid()) received $(length(sc[1,1].cuts)) cuts")
#     m.stagecuts = deepcopy(sc)
#     oldn = Array(Int, (m.stages, m.markov_states,2))
#     for stage=1:m.stages
#         for markovstate=1:m.markov_states
#             oldn[stage, markovstate, 1] = length(m.stagecuts[stage,markovstate].cuts)
#             oldn[stage, markovstate, 2] = length(m.stagecuts[stage,markovstate].samplepoints)
#         end
#     end
#     rebuild_stageproblems!(m)
#     obj = zeros(n)
#     for i=1:n
#         obj[i] = backward_pass!(m, true)
#     end
#     N = getN(m.stagecuts[1,1])
#     res = Array(Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, (m.stages, m.markov_states))
#     for stage=1:m.stages
#         for markovstate=1:m.markov_states
#             res[stage, markovstate] = gettuple(m.stagecuts[stage, markovstate], oldn[stage,markovstate,1], oldn[stage,markovstate,2])
#         end
#     end
#     res, obj
# end
# getN{N}(x::StageCuts{N}) = N
# gettuple{N}(sc::StageCuts{N}, oldcuts::Int, oldsamplepoints::Int) = (sc.samplepoints[(oldsamplepoints+1):end], sc.cuts[(oldcuts+1):end])
#
# function merge_stagecuts!{N}(m::SDDPModel, s::StageCuts{N}, stagecuts::Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}})
#     s.samplepoints = unique(union(s.samplepoints, stagecuts[1]))
#     s.cuts = unique(union(s.cuts, stagecuts[2]))
#     recalculate_dominance!(m.sense, s)
# end
#
# function reduce_async_backwards_pass!{N}(m::SDDPModel, results::Array{Tuple{Vector{NTuple{N,Float64}}, Vector{Cut{N}}}, 2})
#     for stage=1:m.stages
#         for markovstate=1:m.markov_states
#             merge_stagecuts!(m, m.stagecuts[stage, markovstate], results[stage, markovstate])
#         end
#     end
#     m
# end
#
# function update_model!(m::SDDPModel, results, obj)
#     reduce_async_backwards_pass!(m, results)
#     rebuild_stageproblems!(m)
#
#     # Recalculate bound
#     solve_all_stage_problems!(m, 1)
#     set_valid_bound!(m)
#
#     test_and_set_ci!(m, obj)
# end
#
# function async_solve{N1<:Real, N2<:Real}(mm::SDDPModel, simulation_passes::Int, convergence_test_frequency::Int, maximum_iterations::Int, beta_quantile::N1, risk_lambda::N2,  cut_selection_frequency::Int, convergence_test::Bool, cuts_per_processor::Int)
#     # Intialise model on worker processors
#     nworkers = initialise_workers!(mm)
#
#     # Initialise cut selection storage since we need it to communicate between processors
#     initialise_cut_selection!(mm)
#
#     # Set risk aversion parameters
#     mm.beta_quantile, mm.risk_lambda = beta_quantile, risk_lambda
#
#     total_iterations = 0
#     objectives = RollingStorage(convergence_test_frequency)
#
#     function updatemodel!(results, obj)
#         total_iterations+=cuts_per_processor
#
#         # Update simulated bound
#         insert!(objectives, obj)
#
#         # Combine cuts into model
#         update_model!(mm, results, objectives)
#
#         # Output to user
#         # print_stats(m, objectives.n)
#         println("$(total_iterations), $(mm.valid_bound), $(mm.confidence_interval)")
#
#         return
#     end
#     get_stagecuts() = (mm.stagecuts)
#     function check_continue()
#         total_iterations < maximum_iterations && !terminate(rtol(mm) < 0., convergence_test)
#     end
#     @sync begin
#         for procid in workers()
#             @async begin
#                 while check_continue()
#                     stagecuts = get_stagecuts()
#                     results, obj = remotecall_fetch(async_worker_backward_pass!, procid, stagecuts, cuts_per_processor)
#                     updatemodel!(results, obj)
#                 end
#             end
#         end
#     end
#
#     if terminate(rtol(mm) < 0., convergence_test)
#         return CONVERGENCE_TERMINATION
#     else
#         return ITERATION_TERMINATION
#     end
# end
#
# type RollingStorage{T,N}
#     values::Vector{T}   # values
#     i::Int              # next index
#     n::Int              # num values
#     sum::Float64        # mean of stored values
#     sumsquares::Float64 # sum of squares of stored values
# end
# RollingStorage(T::DataType, n::Int) = RollingStorage{T,n}(zeros(T, n), 1, 0, 0., 0.)
# RollingStorage(n::Int) = RollingStorage(Float64, n)
#
# function Base.insert!{T,N}(s::RollingStorage{T,N}, x::T)
#     s.sum += (x - s.values[s.i])            # Update sum
#     s.sumsquares += (x^2- s.values[s.i]^2)  # Update sum of squares
#     s.values[s.i] = x                       # Store value
#     s.i += 1                                # Increment index
#     if s.i > N
#         s.i = 1                      # Wrap around if necessary
#     end
#     if s.n < N
#         s.n += 1                     # Update total values stored
#     end
#     return
# end
# function Base.insert!{T,N}(s::RollingStorage{T,N}, X::Vector{T})
#     for x in X
#         insert!(s, x)
#     end
# end
# Base.mean{T,N}(s::RollingStorage{T,N}) = s.sum / s.n
# function Base.std{T,N}(s::RollingStorage{T,N})
#     μ = mean(s)
#     sqrt((s.sumsquares - 2 * μ * s.sum + s.n * μ^2) / (s.n-1))
# end
# Base.length{T,N}(s::RollingStorage{T,N}) = s.n
# beta_quantile{T,N}(s::RollingStorage{T,N}, beta) = beta_quantile(s.values[1:s.n], beta)
#
# function historicalsimulation(m::SDDPModel, vars::Vector{Symbol}=Symbol[]; markov=ones(Int, m.stages), kwargs...)
#     n=1
#     results = Dict{Symbol, Any}(:Objective=>zeros(Float64,n))
#     for (s, t) in vcat(collect(zip(vars, fill(Any, length(vars)))), [(:Scenario, Int), (:Markov, Int), (:Future, Float64), (:Current, Float64)])
#         results[s] = Array(Vector{t}, m.stages)
#         for i=1:m.stages
#             results[s][i] = Array(t, n)
#         end
#     end
#
#     scenario = 0
#
#     # Initialise the objective storage
#     obj = 0.
#
#     # For all stages
#     for stage=1:m.stages
#         # realise on scenario
#         sp = m.stage_problems[stage, markov[stage]]
#         scenario = load_scenario!(m, sp)
#         data = stagedata(sp)
#         for (key, series) in kwargs
#             @assert haskey(data.scenario_constraint_names, key)
#             cidx = data.scenario_constraint_names[key]
#             JuMP.setRHS(data.scenario_constraints[cidx][1], series[stage])
#         end
#
#         # solve
#         solve!(sp)
#
#         # Add objective (stage profit only)
#         obj += get_true_value(sp)
#
#         # Save results if necesary
#         store_results!(results, vars, sp, stage, 1, markov[stage], scenario)
#         # pass forward if necessary
#         if stage < m.stages
#             pass_states_forward!(m, stage, markov[stage])
#         end
#     end
#
#     results
# end
#
# """
# Calculate Expected objective
# """
# function cvar(x, beta::Float64=1., lambda::Float64=1.)
#     @assert beta >= 0 && beta <= 1.
#     @assert lambda >= 0 && lambda <= 1.
#     cv = lambda * mean(x) + (1 - lambda) * beta_quantile(x, beta)
#     return (cv, cv)
# end
# function beta_quantile(x, beta)
#     mean(x[x.<quantile(x, beta)])
# end
#
# # @variable(m, x[1:10])
# # @variable(m, y)
# # @variable(m, z[1:2, [:a, :b]])
# #
# # results = simulate(m, 1000, x, y, z, z[1,:a] bystage=true)
# # results[x][...key...][stage] >> [ ... realisations ...]
# # results[y][stage] >> [ ... realisations ...]
# # results[z][...key...][stage] >> [ ... realisations ...]
# # results[z[1,:a]][stage] >> [ ... realisations ...]
# #
# # results = simulate(m, 1000, x, y, z, z[1,:a] bystage=false)
# # results[x][...key...][...realisation...] >> [ ... stages ...]
# # results[y][...realisation...] >> [ ... stages ...]
# # results[z][...key...][...realisation...] >> [ ... stages ...]
# # results[z[1,:a]][...realisation...] >> [ ... stages ...]
# #
# # function simulate{M,N,S,T,SEN}(m::SDDPModel{M,N,S,T,SEN}, n::Int, nargs...; bystage=true)
# #     results = initialise_results!(Int(M), n, nargs, bystage)
# #
# #     markov, old_markov = 0, 0
# #
# #     for pass=1:n
# #         if m.init_markov_state==0
# #             markov = transition(m, 1, m.init_markov_state)
# #         else
# #             markov = m.init_markov_state
# #         end
# #         for stage=1:M
# #             old_markov = markov
# #
# #             scenario = rand(1:S)
# #
# #             sp = m.stage_problems[stage, markov]  # corresponding stage problem
# #             load_scenario!(sp, scenario)
# #
# #             solve!(sp)                           # solve
# #
# #             results[:Objective][pass] += get_true_value(sp)         # Add objective
# #
# #             # pass forward if necessary
# #             if stage < M
# #                 pass_states_forward!(m, stage, old_markov)
# #                 # transition to new scenario
# #                 markov = transition(m, stage, old_markov)
# #             end
# #
# #             store_results!(results, sp, stage, pass, nargs, bystage)
# #         end
# #     end
# #     return results
# # end
# # function initialise_results!(M::Int, n::Int, nargs, bystage)
# #     results = Dict{Any, Any}(:Objective=>zeros(Float64,n))
# #
# #     M = Vector{Float64}[]
# #     if bystage
# #         M = Vector{Float64}[zeros(M) for i=1:n]
# #     else
# #         M = Vector{Float64}[zeros(n) for i=1:M]
# #     end
# #
# #     for jump_container in nargs
# #         if isa(jump_container, JuMP.Variable)
# #             results[jump_container] = copy(M)
# #         elseif isa(jump_container, Array{JuMP.Variable})
# #             for i in eachindex(jump_container)
# #                 results[jump_container][ind2sub(size(jump_container), i)...]=copy(M)
# #             end
# #         elseif isa(jump_container, JuMP.JuMPArray)
# #             for key in keys(jump_container)
# #                 results[jump_container][key...] = copy(M)
# #             end
# #         end
# #     end
# #
# #     return results
# # end
# #
# # function store_results!(results, sp, stage, pass, nargs, bystage)
# #     for jump_container in nargs
# #         if isa(jump_container, JuMP.Variable)
# #             if bystage
# #                 results[jump_container][stage][pass] = getvalue(sp, jump_container)
# #             else
# #                 results[jump_container][pass][stage] = getvalue(sp, jump_container)
# #             end
# #         elseif isa(jump_container, Array{JuMP.Variable})
# #             for i in eachindex(jump_container)
# #                 tup = ind2sub(size(jump_container), i)...
# #                 if bystage
# #                     results[jump_container][tup...][stage][pass] = getvalue(sp, jump_container[tup...])
# #                 else
# #                     results[jump_container][tup...][pass][stage] = getvalue(sp, jump_container[tup...])
# #                 end
# #             end
# #         elseif isa(jump_container, JuMP.JuMPArray)
# #             for key in keys(jump_container)
# #                 if bystage
# #                     results[jump_container][key...][stage][pass] = getvalue(sp, jump_container[key...])
# #                 else
# #                     results[jump_container][key...][pass][stage] = getvalue(sp, jump_container[key...])
# #                 end
# #             end
# #         end
# #     end
# # end
