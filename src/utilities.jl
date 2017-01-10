# Copyright 2017, Oscar Dowson

dominates(::Maximisation, x, y) = x < y
dominates(::Minimisation, x, y) = x > y

ext(m::JuMP.Model) = m.ext[:subproblem]::SubproblemExt

_setobjective!(m::JuMP.Model, obj) = JuMP.setobjective(m, JuMP.getobjectivesense(m), obj)

stageproblem(m::SDDPModel, t::Int, i::Int) = m.stageproblems[t][i]
stageproblem(m::SDDPModel, t::Int) = stageproblem(m, t, 1)

getstate(m::JuMP.Model) = map(getvalue, ext(m).states)

cutstorage(m::SDDPModel, t::Int, i::Int) = m.cutstorage[t][i]
cutstorage(m::SDDPModel, t::Int) = cutstorage(m, t, 1)

numstages(m::SDDPModel) = length(m.stageproblems)
nummarkovstates(m::SDDPModel, stage::Int) = length(m.stageproblems[stage])

numstates(m::JuMP.Model) = length(ext(m).states)
numstates(m::SDDPModel, stage::Int, markovstate::Int) = numstates(stageproblem(stage, markovstate))
numstates(m::SDDPModel, stage::Int) = numstates(m, stage, 1)

numpricescenarios(m::JuMP.Model) = length(ext(m).pricescenarios)
numpricescenarios(m::SDDPModel, stage::Int, markovstate::Int) = numpricescenarios(stageproblem(stage, markovstate))
numpricescenarios(m::SDDPModel, stage::Int) = numpriceribs(m, stage, 1)

numpriceribs(m::JuMP.Model) = length(ext(m).theta)
numpriceribs(m::SDDPModel, stage::Int, markovstate::Int) = numpriceribs(stageproblem(stage, markovstate))
numpriceribs(m::SDDPModel, stage::Int) = numpriceribs(m, stage, 1)

numscenarios(m::JuMP.Model) = length(ext(m).scenarios)
numscenarios(m::SDDPModel, stage::Int, markovstate::Int) = numscenarios(stageproblem(stage, markovstate))
numscenarios(m::SDDPModel, stage::Int) = numscenarios(m, stage, 1)

nummarkovstates(T::Array{Float64, 2}, t::Int) = size(T, 1)
nummarkovstates(T::Vector{Array{Float64, 2}}, t::Int) = size(T[t], 1)

getscenariovec(x::Int, t::Int, i::Int) = ones(x) / x
getscenariovec{T<:AbstractVector}(x::T, t::Int, i::Int) = x
getscenariovec{T<:AbstractVector}(x::Vector{T}, t::Int, i::Int) = x[t]
getscenariovec{T<:AbstractVector}(x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]

getel{T}(x::T, t::Int, i::Int) = x
getel{T}(x::Vector{T}, t::Int, i::Int) = x[t]
getel{T}(x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]

transition(x::Array{Float64, 2}, t, i, j) = x[i, j]
transition(x::Vector{Array{Float64, 2}}, t, i, j) = x[t][i, j]
transition(x, t, i, j) = 1.0

function getsense(x::Symbol)
    if x == :Min
        return Minimisation
    elseif x == :Max
        return Maximisation
    else
        error("The sense $(x) must be one of [:Min, :Max]")
    end
end

JuMP.getdual(x::StateVariable) = getdual(x.con)
JuMP.getvalue(x::StateVariable) = getvalue(x.x)
