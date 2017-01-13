# Copyright 2017, Oscar Dowson

typealias LinearConstraint JuMP.ConstraintRef{JuMP.Model, JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64, JuMP.Variable}}}

abstract AbstractCutOracle
abstract AbstractRiskMeasure

typealias Sense Union{Type{Val{:Min}}, Type{Val{:Max}}}
typealias Minimisation Val{:Min}
typealias Maximisation Val{:Max}

immutable Cut
    intercept::Float64
    coefficients::Vector{Float64}
end

immutable CutContainer
    cut::Cut
    stage::Int
    markov::Int
    rib::Int
end

# A state variable with the dummy dual constraint
immutable StateVariable
    v::JuMP.Variable
    c::LinearConstraint
end

# A binding of constraint and possible RHS
immutable ConstraintRHS
    constraint::LinearConstraint
    rhs::Float64
end

# A list of constraints with the RHS that occurs with some probability
immutable Scenario
    probability::Float64
    arr::Vector{ConstraintRHS}
end
Scenario(p) = Scenario(p, ConstraintRHS[])

immutable Rib
    price::Float64
    theta::JuMP.Variable
end

immutable PriceScenario
    probability::Float64
    dynamics::Function
end

type Subproblem
    model::JuMP.Model
    cutoracle::AbstractCutOracle
    riskmeasure::AbstractRiskMeasure
end
model(s::Subproblem) = s.model

immutable SubproblemExtension
    states::Vector{StateVariable}
    scenarios::Vector{Scenario}
    pricescenarios::Vector{PriceScenario}
    ribs::Vector{Rib}
    transitionprobability::Vector{Float64}
end
SubproblemExtension() = SubproblemExtension(
StateVariable[],
Scenario[],
PriceScenario[],
Rib[],
Float64[]
)

type BackwardPassStore
    obj::Float64
    probability::Float64
    pi::Vector{Float64}
end

type ForwardPassStore
    markov::Int
    price::Float64
    state::Vector{Float64}
end

type Stage
    # list of markov subproblems
    arr::Vector{Subproblem}
    # keep track of all states visited in this stage
    statesvisited::Vector{Vector{Float64}}
    # used when solving all problems in a stage (in a single pass)
    backwardpassstore::Vector{BackwardPassStore}
    # store by pass number (i.e. for scenario incrementation)
    forwardpassstore::Vector{ForwardPassStore}
end
Stage() = Stage(Subproblem[], Vector{Float64}[], BackwardPassStore[], ForwardPassStore[])

type SDDPModel
    sense::Sense
    # the problems
    stageproblems::Vector{Stage}
    # the build function
    buildsubproblem!::Function
    # store all cuts
    cutcontainer::Vector{CutContainer}
end
SDDPModel(sense, buildsubproblem!::Function) = SDDPModel(sense, Stage[], buildsubproblem!, CutContainer[])
