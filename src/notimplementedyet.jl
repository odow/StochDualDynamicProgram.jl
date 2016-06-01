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
# @variable(m, x[1:10])
# @variable(m, y)
# @variable(m, z[1:2, [:a, :b]])
#
# results = simulate(m, 1000, x, y, z, z[1,:a] bystage=true)
# results[x][...key...][stage] >> [ ... realisations ...]
# results[y][stage] >> [ ... realisations ...]
# results[z][...key...][stage] >> [ ... realisations ...]
# results[z[1,:a]][stage] >> [ ... realisations ...]
#
# results = simulate(m, 1000, x, y, z, z[1,:a] bystage=false)
# results[x][...key...][...realisation...] >> [ ... stages ...]
# results[y][...realisation...] >> [ ... stages ...]
# results[z][...key...][...realisation...] >> [ ... stages ...]
# results[z[1,:a]][...realisation...] >> [ ... stages ...]
#
# function simulate{M,N,S,T,SEN}(m::SDDPModel{M,N,S,T,SEN}, n::Int, nargs...; bystage=true)
#     results = initialise_results!(Int(M), n, nargs, bystage)
#
#     markov, old_markov = 0, 0
#
#     for pass=1:n
#         if m.init_markov_state==0
#             markov = transition(m, 1, m.init_markov_state)
#         else
#             markov = m.init_markov_state
#         end
#         for stage=1:M
#             old_markov = markov
#
#             scenario = rand(1:S)
#
#             sp = m.stage_problems[stage, markov]  # corresponding stage problem
#             load_scenario!(sp, scenario)
#
#             solve!(sp)                           # solve
#
#             results[:Objective][pass] += get_true_value(sp)         # Add objective
#
#             # pass forward if necessary
#             if stage < M
#                 pass_states_forward!(m, stage, old_markov)
#                 # transition to new scenario
#                 markov = transition(m, stage, old_markov)
#             end
#
#             store_results!(results, sp, stage, pass, nargs, bystage)
#         end
#     end
#     return results
# end
# function initialise_results!(M::Int, n::Int, nargs, bystage)
#     results = Dict{Any, Any}(:Objective=>zeros(Float64,n))
#
#     M = Vector{Float64}[]
#     if bystage
#         M = Vector{Float64}[zeros(M) for i=1:n]
#     else
#         M = Vector{Float64}[zeros(n) for i=1:M]
#     end
#
#     for jump_container in nargs
#         if isa(jump_container, JuMP.Variable)
#             results[jump_container] = copy(M)
#         elseif isa(jump_container, Array{JuMP.Variable})
#             for i in eachindex(jump_container)
#                 results[jump_container][ind2sub(size(jump_container), i)...]=copy(M)
#             end
#         elseif isa(jump_container, JuMP.JuMPArray)
#             for key in keys(jump_container)
#                 results[jump_container][key...] = copy(M)
#             end
#         end
#     end
#
#     return results
# end
#
# function store_results!(results, sp, stage, pass, nargs, bystage)
#     for jump_container in nargs
#         if isa(jump_container, JuMP.Variable)
#             if bystage
#                 results[jump_container][stage][pass] = getvalue(sp, jump_container)
#             else
#                 results[jump_container][pass][stage] = getvalue(sp, jump_container)
#             end
#         elseif isa(jump_container, Array{JuMP.Variable})
#             for i in eachindex(jump_container)
#                 tup = ind2sub(size(jump_container), i)...
#                 if bystage
#                     results[jump_container][tup...][stage][pass] = getvalue(sp, jump_container[tup...])
#                 else
#                     results[jump_container][tup...][pass][stage] = getvalue(sp, jump_container[tup...])
#                 end
#             end
#         elseif isa(jump_container, JuMP.JuMPArray)
#             for key in keys(jump_container)
#                 if bystage
#                     results[jump_container][key...][stage][pass] = getvalue(sp, jump_container[key...])
#                 else
#                     results[jump_container][key...][pass][stage] = getvalue(sp, jump_container[key...])
#                 end
#             end
#         end
#     end
# end
