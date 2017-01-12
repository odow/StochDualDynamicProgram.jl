# Copyright 2017, Oscar Dowson

const comparison_symbols = [:(<=), :(>=), :(==)]
is_comparison(x) = Base.Meta.isexpr(x, :comparison) || (Base.Meta.isexpr(x, :call) && x.args[1] in comparison_symbols)

_setobjective!(m::JuMP.Model, obj) = JuMP.setobjective(m, JuMP.getobjectivesense(m), obj)

#
#   Model.jl
#

model(x::SubproblemExtension) = x.model
states(x::SubproblemExtension) = x.states
scenarios(x::SubproblemExtension) = x.scenarios
pricescenarios(x::SubproblemExtension) = x.pricescenarios
ribs(x::JuMP.Model) = ribs(ext(x))
model(x::JuMP.Model) = model(ext(x))
states(x::JuMP.Model) = states(ext(x))
scenarios(x::JuMP.Model) = scenarios(ext(x))
pricescenarios(x::JuMP.Model) = pricescenarios(ext(x))
ribs(x::JuMP.Model) = ribs(ext(x))


#
#    Set markov transition
#

numberofmarkovstates(x::Int, t::Int) = x
numberofmarkovstates(x::Vector{Int}, t::Int) = x[t]

function settransition!(ex::SubproblemExtension, transition, markovstates, T, t, i)
    if t < T
        J = numberofmarkovstates(markovstates, t+1)
        tmat = extracttransition(transition, t, i, J)
        settransition!(sp, tmat)
    end
    nothing
end

function settransition!(ex::SubproblemExtension, T::Vector{Float64})
    while length(ex.transitionprobability) > 0
        pop!(ex.transitionprobability)
    end
    for t in T
        push!(ex.transitionprobability, t)
    end
end
settransition!(sp::JuMP.Model, T::Vector{Float64}) = settransition!(ext(sp), T)

extracttransition(::Void, t, i, J::Int) = ones(Float64, J) / J
extracttransition(T::Array{Float64, 2}, t, i, J::Int) = T[i, :]
extracttransition(T::Vector{Array{Float64, 2}}, t, i, J::Int) = extracttransition(T[t], t, i, J)

#
#    Set scenarios
#
scenarioprob(x::Void, t, i) = Float64[]
scenarioprob(x::Int, t, i) = ones(Float64, x) / x
scenarioprob(x::Vector{Float64}, t, i) = x
scenarioprob{T<:Union{Int, Vector{Float64}}}(x::Vector{T}, t, i) = scenarioprob(x[t], t, i)
scenarioprob{T<:Union{Int, Vector{Float64}}}(x::Vector{Vector{T}}, t, i) = scenarioprob(x[t][i], t, i)
function setscenarios!(sp::JuMP.Model, scenarios, t, i)
    for s in scenarioprob(scenarios)
        push!(scenarios(sp), Scenario(s))
    end
end
