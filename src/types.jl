# ==============================================================================
#   Optimisation Sense
abstract AbstractSense
immutable Max <: AbstractSense end
immutable Min <: AbstractSense end
typealias Sense Union{Type{Max}, Type{Min}}

# ==============================================================================
#   Cuts
"""
    Cut{N}

A single cut benders cut for N state variables. Of the form

    theta >= intercept + dot(coefficients, x)

Parameters:

    N                the number of state variables (length of coefficients)

Fields:

    intercept        intercept of cut
    coefficients     vector of coefficients
"""
type Cut{N}
    intercept::Float64
    coefficients::Vector{Float64}
end
Cut(intercept, coefficients::Vector) = Cut{length(coefficients)}(intercept, coefficients)
Cut(N) = Cut{N}(0., zeros(N))

"""
    StageCuts{N}

A type to hold all the cuts in a single subproblem with dominance.

Parameters:

    N                  the number of state variables (length of coefficients)

Fields:

    n                  the number of samplepoints
    cuts               a vector of all the discovered cuts
    samplepoints       the sample points visited by the algorithm
    nondominated       the number of samplepoints that the cut i is nondominated at
    activecut          the index of the cut active at samplepoint i
"""
type StageCuts{N}
    n::Int                                    # Number of sample points in stage problem
    samplepoints::Vector{NTuple{N, Float64}}  # x to evaluate at
    cuts::Vector{Cut{N}}                      # list of cuts in stage problem.
    nondominated::Vector{Int}                 # number of points cut is nondominated at length(nondomindated) == length(cuts)
    activecut::Vector{Int}                    # index of cut that is active at point x(i) length(active_cut) == n
end
function StageCuts(sp::Model, bound)
    N = length(stagedata(sp).state_vars)
    StageCuts(0, NTuple{N, Float64}[], Cut{N}[Cut(bound, zeros(N))], Int[0], Int[])
end

"""
    StageDataExt

JuMP extension structure for SDDPModel subproblems.

Fields:

    state_vars                 Vector of state variables in the subproblem
    dual_constraints           Vector of dummy constraints to get incoming duals
    theta                      Cost/Value to go variable

    last_markov                Index of the last markov state visited

    scenario_constraints       Vector of scenario constraints
    scenario_constraint_names  Dictionary where key is symbol of scenarioconstaint name and value is integer index of constraint in scenario constraints

    objective_values           Vector of objective values for each scenario
    dual_values                Vector of vectors of dual values by scenario

    beta_quantile              Beta quantile for risk aversion
    lambda_weight              Weighing on expectation for risk aversion
    weightings_matrix          Matrix for cut weights

    stage_profit               Expression for the stage profit

    regularisecoefficient
    regularisecons
    regularisepen
"""
type StageDataExt
    # vector of state variables
    state_vars::Vector{JuMP.Variable}
    # vector of dummy constraints to get duals from
    dual_constraints::Vector{JuMP.ConstraintRef}
    # value/cost-to-go variable
    theta::Union{Void, JuMP.Variable}

    # index of last markov state
    last_markov::Int

    # Storage for scenarioconstraints
    scenario_constraints::Vector{Tuple{Any, Vector{Any}}}
    scenario_constraint_names::Dict{Symbol, Int}

    # Store values across scenarios
    objective_values::Vector{Float64}
    dual_values::Vector{Vector{Float64}}

    beta_quantile::Float64
    lambda_weight::Float64
    weightings_matrix::Array{Float64, 2}

    stage_profit

    # Regularisation stuff
    regularisecoefficient::Float64
    regularisecons::Vector{JuMP.ConstraintRef}
    regularisepen::Union{Void, JuMP.Variable}
end
StageDataExt(scenarios::Int=1) = StageDataExt(
    Variable[],
    ConstraintRef[],
    nothing,

    0,

    Tuple{Any, Vector{Any}}[],
    Dict{Symbol, Int}(),

    zeros(scenarios),
    Vector{Float64}[],

    1.,
    1.,
    Array(Float64, (0,0)),

    0.,

    1.,
    ConstraintRef[],
    nothing
    )

# ==============================================================================
#   Forward Pass
"""
    ForwardPassData

Storage for the forward pass of the algorithm

Fields:

    n       number of realisations in forward pass
    x       array of states visited x[state, stage, realisation]
    obj     vector of objective values obj[realisation]
    W       array of Markov states visited W[stage, realisation]
"""
type ForwardPassData
    n::Int
    x::Vector{Array{Float64, 2}}
    obj::Vector{Float64}
    W::Vector{Vector{Int}}
end
ForwardPassData() = ForwardPassData(
                    0,
                    Array{Float64, 2}[],
                    Float64[],
                    Vector{Int}[]
                    )

"""
    SDDPModel{T, M, S, X, TM}

The SDDP model object.

Parameters:

    T                          The number of stages in the problems
    M                          The number of markov states in the model
    S                          The number of scenarios in the model
    X                          The sense (max/min) of the model
    TM                         The type of markov transition matrix

Fields:

    stage_problems             Array of JuMP models for each subproblem (stage x markov state)
    transition                 The transition matrix
    initial_markov_state       The initial markov state
    scenario_probability       Scenario probability support vector
    confidence_interval        Estimated value of policy
    valid_bound                Valid bound for value of policy
    build_function!            Subproblem generation function
    stagecuts                  Array for the stage cuts
    solver                     MathProgBase solver for subproblems
    valuetogobound             Bound on the cost/value to go
    forwardstorage             Storage for the forward pass
"""
type SDDPModel{T, M, S, X<:AbstractSense, TM}
    # array of JuMP subproblems
    stage_problems::Array{JuMP.Model, 2}

    # markov transition matrix
    transition::Array{Float64, TM}
    # initial markov state
    initial_markov_state::Int

    # scenario probability support vector
    scenario_probability::WeightVec

    # approximate statistical bound for problem
    confidence_interval::Tuple{Float64, Float64}
    valid_bound::Float64

    # function that builds the stage problems
    build_function!::Function

    # Array for all the stage cuts
    stagecuts::Array{StageCuts, 2}

    solver::MathProgBase.AbstractMathProgSolver
    valuetogobound::Float64

    forwardstorage::ForwardPassData
end

# Define some short cuts to access members
stagecut(m::SDDPModel, t, i) = m.stagecuts[t,i]
subproblem(m::SDDPModel, t, i) = m.stage_problems[t,i]
subproblems(m::SDDPModel, t) = m.stage_problems[t,:]
stagedata(m::SDDPModel, t, i) = stagedata(subproblem(m,t,i))

# ==============================================================================
#   Risk Measures
abstract RiskMeasure

immutable NestedCVar <: RiskMeasure
    beta::Float64
    lambda::Float64
    NestedCVar{T<:Real, S<:Real}(beta::T, lambda::S) = (@assert beta >= 0. && beta <= 1. && lambda >= 0. && lambda <= 1.; new(beta, lambda))
end

NestedCVar(;beta=1., lambda=1.) = NestedCVar(beta, lambda)
Expectation() = NestedCVar(1., 1.)

# ==============================================================================
#   Regularisation
abstract Regularisation

immutable NoRegularisation <: Regularisation
    initial::Float64
    decayrate::Float64
    NoRegularisation() = new(0., 0.)
end

immutable LinearRegularisation <: Regularisation
    initial::Float64
    decayrate::Float64
    LinearRegularisation(initial=1., decayrate=0.95) = new(initial, decayrate)
end

immutable QuadraticRegularisation <: Regularisation
    initial::Float64
    decayrate::Float64
    QuadraticRegularisation(initial=1., decayrate=0.95) = new(initial, decayrate)
end

# ==============================================================================
#   Cut Selection
abstract CutSelectionMethod

immutable LevelOne <: CutSelectionMethod
    frequency::Int
    LevelOne(frequency::Int) = (@assert frequency > 0; new(frequency))
end

immutable Deterministic <: CutSelectionMethod
    frequency::Int
    Deterministic(frequency::Int) = (@assert frequency > 0; new(frequency))
end

immutable NoSelection <: CutSelectionMethod
    frequency::Int
    NoSelection() = new(0)
end

# ==============================================================================
#   Parallel
abstract AbstractParallel
immutable Serial <: AbstractParallel end
immutable ConvergenceTest <: AbstractParallel end
immutable BackwardPass <: AbstractParallel
    cuts_per_processor::Int
    BackwardPass(n::Int) = (@assert n > 0; new(n))
end
BackwardPass() = BackwardPass(1)

type Parallel{AP1<:Union{ConvergenceTest, Serial}, AP2<:Union{BackwardPass, Serial}}
    convergence_test::AP1
    backward_pass::AP2
end
Parallel(bp::BackwardPass, fp::ConvergenceTest) = Parallel(fp, bp)
Parallel(bp::BackwardPass) = Parallel(Serial(), bp)
Parallel(fp::ConvergenceTest) = Parallel(fp, Serial())
Parallel() = Parallel(Serial(), Serial())

# ==============================================================================
#   Terminsation status
const BOUND_TERMINATION     = :BoundConvergence
const POLICY_TERMINATION    = :PolicyConverence
const ITERATION_TERMINATION = :MaximumIterations
const UNKNOWN_TERMINATION   = :Unknown

# ==============================================================================
#   Solution
type SolutionLog
    ci_lower::Float64
    ci_upper::Float64
    bound::Float64
    cuts::Int
    time_backwards::Float64
    simulations::Int
    time_forwards::Float64
    time_cutselection::Float64
end
SolutionLog() = SolutionLog(0., 0., 0., 0, 0., 0, 0., 0.)
Base.copy(s::SolutionLog) = SolutionLog(s.ci_lower, s.ci_upper, s.bound, s.cuts, s.time_backwards, s.simulations, s.time_forwards, s.time_cutselection)

type Solution
    status::Symbol
    iterations::Int
    trace::Vector{SolutionLog}
end

Solution() = Solution(UNKNOWN_TERMINATION, 0, SolutionLog[])

setStatus!(s::Solution, sym::Symbol) = (s.status = sym)
