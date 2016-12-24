#  Copyright 2016, Oscar Dowson

isdefined(Base, :__precompile__) && __precompile__()

module StochDualDynamicProgram

using JuMP
using MathProgBase, Clp
using Distributions, StatsBase
using JSON
using Compat

import Base.dot

export SDDPModel,
    # macros
    @state,@states,
    @scenarioconstraint, @scenarioconstraints,
    @relaxedconstraint, @relaxedconstraints,
    @stageprofit,
    @visualise,
    # Model functions
    simulate, historicalsimulation, loadcuts!,
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
    # Backward Pass options
    BackwardPass,
    # Regularisation options
    NoRegularisation, LinearRegularisation, QuadraticRegularisation,
    # query solve attributes
    status,
    # non convex helpers
    SOSII!

include("types.jl")
include("macros.jl")
include("model.jl")
include("forwardpass.jl")
include("risk_aversion.jl")
include("backwardpass.jl")
include("cut_selection.jl")
include("simulate.jl")
include("visualiser/visualise.jl")
include("parallel.jl")
include("MIT_licencedcode.jl")
include("print.jl")
include("expectedvalueproblem.jl")
include("nonconvex.jl")

const PRINTALL   = 4
const PRINTINFO  = 3
const PRINTTRACE = 2
const PRINTTERM  = 1
const PRINTNONE  = 0

"""
    solve(m[; kwargs...])

Solve the SDDPModel
"""
function JuMP.solve{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM};
    maximum_iterations      = 1,
    expectedvalueiterations = 0,
    policy_estimation       = MonteCarloEstimator(),
    bound_convergence       = BoundConvergence(),
    forward_pass            = ForwardPass(),
    backward_pass           = BackwardPass(),
    risk_measure            = Expectation(),
    cut_selection           = NoSelection(),
    parallel                = Serial(),
    output                  = nothing,
    cut_output_file         = nothing,
    print_level             = PRINTALL
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

    print_level >= PRINTTRACE && printheader()

    resizeforwardstorage!(m, 1)
    if expectedvalueiterations > 0
        print_level >= PRINTINFO && info("Running Expected Value Iterations")
        if risk_measure.lambda < (1-1e-5)
            warn("Cannot run expected value problems iterations with risk aversion")
        else
            for evit = 1:expectedvalueiterations
                solution.iterations += 1
                expectedforwardpass!(log, m, 1, cut_selection, forward_pass)
                setconfidenceinterval!(log, m, 0.95)
                expectedbackwardpass!(log, m, risk_measure, forward_pass.regularisation, backward_pass)
                notconverged = setbound!(log, m, bound_convergence)

                # print solution to user
                print_level >= PRINTTRACE && print(m, log, output, false)
                push!(solution.trace, copy(log))
            end
        end
        print_level >= PRINTINFO && info("Proceeding to normal iterations.")
    end

    notconverged = true
    while notconverged
        # Starting a new iteration
        solution.iterations += 1

        # forward pass
        nscenarios = forwardpass!(log, m, solution.iterations, cut_selection, forward_pass, parallel.forward)

        # estimate bound
        notconverged, ismontecarlo = estimatebound!(log, m, policy_estimation, solution.iterations, parallel.montecarlo, print_level)
        !notconverged && setStatus!(solution, POLICY_TERMINATION)

        if notconverged
            # backward pass
            backwardpass!(log, m, nscenarios, risk_measure, forward_pass.regularisation, parallel.backward, backward_pass)

            # run cut selection
            if parallel.backward
                # We have probably dragged various cuts through here...
                recalculate_dominance!(cut_selection, m)
            end
            cutselection!(log, m, cut_selection, solution.iterations, print_level)

            if cut_output_file != nothing
                # Write cuts to file if appropriate
                writecuts!(m, cut_output_file)
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
        print_level >= PRINTTRACE && print(m, log, output, ismontecarlo)

        push!(solution.trace, copy(log))
    end
    print_level >= PRINTTERM && info("Terminating with status $(solution.status).")
    return solution
end

end
