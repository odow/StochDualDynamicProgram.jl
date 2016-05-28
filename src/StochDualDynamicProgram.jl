isdefined(Base, :__precompile__) && __precompile__()

module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp
using Formatting
using StatsBase

import Base.dot

export SDDPModel,
    @state, @scenarioconstraint, @scenarioconstraints, @stageprofit,
    @visualise,
    simulate, loadcuts!,
    LevelOne, Deterministic, NoSelection,
    ConvergenceTest, BackwardPass, Parallel,
    NestedCVar,
    Convergence,
    NoRegularisation, LinearRegularisation, QuadraticRegularisation

include("types.jl")
include("macros.jl")
include("model.jl")
include("forwardpass.jl")
include("risk_aversion.jl")
include("backwardpass.jl")
include("cut_selection.jl")
include("print.jl")
include("simulate.jl")


# ==============================================================================
#   Convergence
type Convergence
    n
    frequency::Int
    terminate::Bool
    quantile::Float64
    variancereduction::Bool
    function Convergence(simulations, frequency::Int, terminate::Bool=false, quantile::Float64=0.95, variancereduction::Bool=true)
        @assert quantile >= 0 && quantile <= 1.
        new(simulations, frequency, terminate, quantile, variancereduction)
    end
end

Convergence(;simulations=1, frequency=1, terminate=false, quantile=0.95, variancereduction=true) = Convergence(simulations, frequency, terminate, quantile, variancereduction)

"""
    MonteCarloEstimator

Solve option to control behaviour of out-of-sample monte carlo estimates for the policy.

Fields:

    frequency          frequency (by iteration) /w which MonteCarlo estimate is run
    minsamples         minimium number of samples before testing for convergence
    maxsamples         maximum number of samples before testing for convergence
    step               step size for simulation incrementation
    terminate          true = terminate if confidence level contains lower bound
    confidencelevel    size of confidence level to construct to check for convergence
"""
type MonteCarloEstimator
    frequency::Int
    minsamples::Int
    maxsamples::Int
    step::Int
    terminate::Bool
    confidencelevel::Float64
end
MonteCarloEstimator(;
    frequency       = 0,
    minsamples      = 10,
    maxsamples      = minsamples,
    step            = 0,
    terminate       = false,
    confidencelevel = 0.95) = MonteCarloEstimator(
                                frequency,
                                minsamples,
                                maxsamples,
                                step,
                                terminate,
                                confidencelevel
                                )

"""

Fields:

    delta       change in lower bound between iterations
    n           number of iterations with change in lower bound less than delta before terminating
"""
type Bound
    delta::Float64
    n::Int
end

"""
    solve(m[; kwargs...])

Solve the SDDPModel
"""
function JuMP.solve{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM};
    maximum_iterations = 1,
    convergence        = Convergence(1, 1, false, 0.95),
    forward_scenarios  = 1,
    risk_measure       = Expectation(),
    cut_selection      = NoSelection(),
    parallel           = Parallel(),
    regularisation     = NoRegularisation(),
    output             = nothing,
    cut_output_file    = nothing
    )

    setriskmeasure!(m, risk_measure)

    solution = Solution()
    log = SolutionLog()

    printheader()

    cutswrittentofile = 0
    while solution.iterations < maximum_iterations
        # forward pass
        forwardpass!(log, m, forward_scenarios, isa(cut_selection, LevelOne))

        # estimate bound
        setconfidenceinterval!(log, m, 0.95)

        # backward pass
        backwardpass!(log, m, forward_scenarios, risk_measure, regularisation)

        # Calculate a new upper bound
        setbound!(log, m)

        # run cut selection
        cutselection!(log, m, cut_selection, solution.iterations)

        if cut_output_file != nothing
            # Write cuts to file if appropriate
            writecuts!(m, cut_output_file, cutswrittentofile)
            cutswrittentofile += forward_scenarios
        end

        # print solution to user
        print(m, log, output)

        push!(solution.trace, copy(log))
        solution.iterations += 1
    end
    return solution
end

# a wrapper for forward pass timings
function forwardpass!(log::SolutionLog, m::SDDPModel, n, storesamplepoints)
    tic()
    forwardpass!(m, n, storesamplepoints)
    log.simulations += n
    log.time_forwards += toq()
    return
end

# a wrapper for backward pass timings
function backwardpass!(log::SolutionLog, m::SDDPModel, forward_scenarios, risk_measure, regularisation)
    tic()
    backwardpass!(m, risk_measure, regularisation)
    log.cuts += forward_scenarios
    log.time_backwards += toq()
    return
end

function setbound!(log::SolutionLog, m::SDDPModel)
    setbound!(m)
    log.bound = m.valid_bound
    return
end

# a wrapper for cut selection timings
function cutselection!(log::SolutionLog, m::SDDPModel, cutselection, iterations)
    tic()
    cutselection!(m, cutselection, iterations)
    log.time_cutselection += toq()
    return
end

# a wrapper for estimating the confidence interval
function setconfidenceinterval!(log, m, conflevel)
    setconfidenceinterval!(m, estimatebound(getobj(m), conflevel))
    log.ci_lower, log.ci_upper = m.confidence_interval
    return
end

end
