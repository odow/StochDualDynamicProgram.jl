# Copyright 2017, Oscar Dowson

function checkzerotoone(x)
    @assert x <= 1
    @assert x >= 0
end

dominates(::Maximisation, x, y) = x < y
dominates(::Minimisation, x, y) = x > y

ext(m::JuMP.Model) = m.ext[:subproblem]::SubproblemExt

_setobjective!(m::JuMP.Model, obj) = JuMP.setobjective(m, JuMP.getobjectivesense(m), obj)

stageproblem(m::SDDPModel, t::Int, i::Int) = m.stageproblems[t][i]
stageproblem(m::SDDPModel, t::Int) = stageproblem(m, t, 1)

function getstate(m::JuMP.Model)
    y = zeros(numstates(m))
    getstate!(m, y)
    y
 end
function getstate!(m::JuMP.Model, x::Vector{Float64})
    @assert length(x) == numstates(m)
    for i=1:numstates(m)
        x[i] = getvalue(state(m, i))
    end
end
state(m::JuMP.Model, i::Int) = ext(m).states[i]

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
transition(x::Vector{Array{Float64, 2}}, t, i, j) = x[t+1][i, j]
transition(x, t, i, j) = 1.0
transition(m::SDDPModel, t, i, j) = transition(m.transition, t, i, j)
function transition(m::SDDPModel, t, i)
    r = rand()
    for j = 1:nummarkovstates(m, t)
        r -= transition(m, t, i, j)
        if r <= 0.0
            return j
        end
    end
    return -1
end

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

addstatevisited!(m, t, x) = push!(m.statesvisited[t], copy(x))
probability(scenario::Scenario) = scenario.probability
probability{T}(price::PriceScenario{T}) = price.probability

storeobj!(m, idx, obj) = (m.backwardstorage[idx].obj = obj)
storeprob!(m, idx, prob) = (m.backwardstorage[idx].probability = prob)
dualstore(m, idx) = m.backwardstorage[idx].pi

pricescenario(m::JuMP.Model, p) = ext(m).pricescenarios[p]
scenario(m::JuMP.Model, s) = ext(m).scenarios[s]
rib(m::JuMP.Model, r) = ext(m).theta[r]

ForwardStorage(m::SDDPModel) = ForwardStorage(numstages(m), numstates(stageproblem(m, 1, 1)))

getmarkov(m::SDDPModel, pass::Int, stage::Int) = m.forward_storage[pass][stage].markov
getprice(m::SDDPModel, pass::Int, stage::Int) = m.forward_storage[pass][stage].price
getstate(m::SDDPModel, pass::Int, stage::Int) = m.forward_storage[pass][stage].state

setmarkov!(m::SDDPModel, pass::Int, stage::Int, x) = (m.forward_storage[pass][stage].markov = x)
setprice!(m::SDDPModel, pass::Int, stage::Int, x) = (m.forward_storage[pass][stage].price = x)
function setstate!(m::SDDPModel, pass::Int, stage::Int, x)
    for i=1:length(x)
        m.forward_storage[pass][stage].state[i] = x[i]
    end
end
