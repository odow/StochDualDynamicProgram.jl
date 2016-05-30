isdefined(Base, :__precompile__) && __precompile__()

module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp
using Formatting
using Distributions, StatsBase
using JSON

import Base.dot

export SDDPModel,
    # macros
    @state, @scenarioconstraint, @scenarioconstraints, @stageprofit, @visualise,
    # Model functions
    simulate, loadcuts!,
    # Cut selection options
    LevelOne, Deterministic, NoSelection,
    # Parallel options
    Serial, Parallel,
    # Risk averse options
    NestedCVar, Expectation,
    # Termination options
    MonteCarloEstimator, BoundConvergence,
    # Forward Pass options
    ForwardPass,
    # Regularisation options
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
include("visualiser/visualise.jl")
include("parallel.jl")

"""
    solve(m[; kwargs...])

Solve the SDDPModel
"""
function JuMP.solve{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM};
    maximum_iterations = 1,
    convergence        = MonteCarloEstimator(),
    bound_convergence  = BoundConvergence(),
    forward_pass       = ForwardPass(),
    risk_measure       = Expectation(),
    cut_selection      = NoSelection(),
    parallel           = Serial(),
    output             = nothing,
    cut_output_file    = nothing
    )

    setriskmeasure!(m, risk_measure)

    if isa(parallel, Parallel)
        if length(workers()) < 2
            warn("Paralleisation requested but Julia is only running with a single worker. Start julia with `julia -p N` or use the `addprocs(N)` function to load N workers. Running in serial mode.")
            parallel = Serial()
        else
            nworkers = initialise_workers!(m)
        end
    end

    solution = Solution()
    log = SolutionLog()

    printheader()

    cutswrittentofile = 0
    notconverged = true
    while notconverged
        # Starting a new iteration
        solution.iterations += 1

        # number of scenarios to sample on forward pass
        nscenarios = getscenarios(forward_pass.scenarios, solution.iterations)

        # forward pass
        forwardpass!(log, m, nscenarios, cut_selection, forward_pass)

        # estimate bound
        notconverged, ismontecarlo = estimatebound!(log, m, convergence, solution.iterations, parallel.montecarlo)
        !notconverged && setStatus!(solution, POLICY_TERMINATION)

        if notconverged
            # backward pass
            backwardpass!(log, m, nscenarios, risk_measure, forward_pass.regularisation)

            # run cut selection
            cutselection!(log, m, cut_selection, solution.iterations)

            if cut_output_file != nothing
                # Write cuts to file if appropriate
                writecuts!(m, cut_output_file, cutswrittentofile)
                cutswrittentofile += nscenarios
            end

            # Calculate a new upper bound
            notconverged = setbound!(log, m, bound_convergence)
            !notconverged && setStatus!(solution, BOUND_TERMINATION)

            if solution.iterations == maximum_iterations
                notconverged = false
                setStatus!(solution, ITERATION_TERMINATION)
            end
        end

        # print solution to user
        print(m, log, output, ismontecarlo)

        push!(solution.trace, copy(log))
    end
    info("Terminating with status $(solution.status).")
    return solution
end

"""
    getscenarios(forwardscenarios, iteration)

This function gets the number of scenarios to sample in iteration `iteration`.
"""
getscenarios(forwardscenarios::Int, iteration::Int) = forwardscenarios
function getscenarios(forwardscenarios::AbstractArray{Int,1}, iteration::Int)
    if iteration < length(forwardscenarios)
        return forwardscenarios[iteration]
    else
        return forwardscenarios[end]
    end
end

function setbound!(log::SolutionLog, m::SDDPModel, bound_convergence)
    old_bound = m.valid_bound
    setbound!(log, m)
    if bound_convergence.after > 0
        if abs(old_bound - m.valid_bound) < bound_convergence.tol
            bound_convergence.n += 1
            if bound_convergence.n > bound_convergence.after
                return false # converged
            end
        else
            bound_convergence.n = 0
        end
    end
    return true # notcongerged
end

"""
    estimatebound!(log, SDDPmodel, convergence, iteration)

Estimate the value of the policy by monte-carlo simulation and using sequential sampling.
"""
function estimatebound!(log::SolutionLog, m::SDDPModel, convergence, iteration::Int, isparallel::Bool)
    if isparallel
        estimatebound!(log, m, convergence, iteration, parallelmontecarloestimation)
    else
        estimatebound!(log, m, convergence, iteration, montecarloestimation)
    end
end

function estimatebound!(log::SolutionLog, m::SDDPModel, convergence, iteration, montecarlofunction::Function)
    notconverged = true
    ismontecarlo = false
    if convergence.frequency > 0 && mod(iteration, convergence.frequency) == 0
            tic()
            info("Running out-of-sample Monte Carlo simulation")
            ismontecarlo = true
            obj = copy(getobj(m))
            if length(obj) < convergence.minsamples
                log.simulations += convergence.minsamples - length(obj)
                push!(obj,
                    montecarlofunction(convergence.antithetic, m,
                        convergence.minsamples - length(obj)
                        )...
                )
            end
            ci = estimatebound(obj, convergence.confidencelevel)
            while isconverged(m, ci, m.valid_bound)
                if length(obj) >= convergence.maxsamples
                    if convergence.terminate
                        notconverged = false
                    end
                    break
                end
                oldlength = length(obj)
                push!(obj,
                    montecarlofunction(convergence.antithetic, m,
                        min(
                            convergence.maxsamples-length(obj),
                            convergence.step
                            )
                        )...
                    )
                log.simulations += length(obj) - oldlength
                ci = estimatebound(obj, convergence.confidencelevel)
            end
            setconfidenceinterval!(m, ci)
            log.ci_lower, log.ci_upper = m.confidence_interval
            log.time_forwards += toq()
    else
        setconfidenceinterval!(log, m, 0.95)
    end
    return notconverged, ismontecarlo
end

# a wrapper for forward pass timings
function forwardpass!(log::SolutionLog, m::SDDPModel, n, cutselection, regularisation)
    tic()
    forwardpass!(m, n, cutselection, regularisation)
    log.simulations += n
    log.time_forwards += toq()
    return
end

# a wrapper for backward pass timings
function backwardpass!(log::SolutionLog, m::SDDPModel, nscenarios, risk_measure, regularisation)
    tic()
    backwardpass!(m, risk_measure, regularisation)
    log.cuts += nscenarios
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
    if cutselection.frequency > 0
        tic()
        cutselection!(m, cutselection, iterations)
        log.time_cutselection += toq()
    end
    return
end

# a wrapper for estimating the confidence interval
function setconfidenceinterval!(log, m, conflevel)
    setconfidenceinterval!(m, estimatebound(getobj(m), conflevel))
    log.ci_lower = log.ci_upper = mean(m.confidence_interval)
    return
end

isconverged(::Type{Min}, ci::NTuple{2, Float64}, bound::Float64) = ci[1] < bound
isconverged(::Type{Max}, ci::NTuple{2, Float64}, bound::Float64) = ci[2] > bound
isconverged{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, ci, bound) = isconverged(X, ci, bound)

end
