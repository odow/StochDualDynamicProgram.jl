# Copyright 2017, Oscar Dowson

const comparison_symbols = [:(<=), :(>=), :(==)]
is_comparison(x) = Base.Meta.isexpr(x, :comparison) || (Base.Meta.isexpr(x, :call) && x.args[1] in comparison_symbols)

_setobjective!(m::JuMP.Model, obj) = JuMP.setobjective(m, JuMP.getobjectivesense(m), obj)

#
#   Model.jl
#

function getsense(x::Symbol)
    if x == :Min
        Minimisation
    elseif x == :Max
        Maximisation
    else
        error("Sense $(x) not supported.")
    end
end

getmodel(x::SubproblemExtension) = x.model
getstates(x::SubproblemExtension) = x.states
getscenarios(x::SubproblemExtension) = x.scenarios
getpricescenarios(x::SubproblemExtension) = x.pricescenarios
getribs(x::SubproblemExtension) = x.ribs
getmodel(x::JuMP.Model) = getmodel(ext(x))
getstates(x::JuMP.Model) = getstates(ext(x))
getscenarios(x::JuMP.Model) = getscenarios(ext(x))
getpricescenarios(x::JuMP.Model) = getpricescenarios(ext(x))
getribs(x::JuMP.Model) = getribs(ext(x))


#
#    Set markov transition
#

numberofmarkovstates(x::Int, t::Int) = x
numberofmarkovstates(x::Vector{Int}, t::Int) = x[t]

function settransition!(ex::SubproblemExtension, transition, markovstates, T, t, i)
    if t < T
        J = numberofmarkovstates(markovstates, t+1)
        tmat = extracttransition(transition, t, i, J)
        settransition!(ex, tmat)
    end
    nothing
end
settransition!(sp::JuMP.Model, transition, markovstates, T, t, i) = settransition!(ext(sp), transition, markovstates, T, t, i)
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
function setscenarios!(sp::JuMP.Model, sc, t, i)
    for s in scenarioprob(sc, t, i)
        push!(getscenarios(sp), Scenario(s))
    end
end
