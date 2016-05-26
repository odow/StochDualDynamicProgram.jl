isdefined(Base, :__precompile__) && __precompile__()

module StochDualDynamicProgram

importall JuMP
using MathProgBase, Clp
using Formatting
using Distributions, StatsBase

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
include("SDDPModel.jl")
include("cut_selection.jl")
include("SDDPalgorithm.jl")
include("parallel.jl")
include("visualiser/visualise.jl")

"""
Solve the model using the SDDP algorithm.
"""
function JuMP.solve(m::SDDPModel; maximum_iterations=1, convergence=Convergence(1, 1, false, 0.95), risk_measure=Expectation(), cut_selection=NoSelection(), parallel=Parallel(), regularisation=NoRegularisation(), output=nothing)

    @assert maximum_iterations >= 0

    if !isa(parallel, Parallel{Serial, Serial}) && length(workers()) < 2
        warn("Paralleisation requested but Julia is only running with a single processor. Start julia with `julia -p N` or use the `addprocs(N)` function to load N workers. Running in serial mode.")
        parallel = Parallel()
    end

    for sp in m.stage_problems
        stagedata(sp).regularisecoefficient = regularisation.initial
    end

    # Initial output for user
    print_stats_header()

    solution = Solution()

    m.QUANTILE = convergence.quantile

    # try
    solve!(m, solution, convergence, maximum_iterations, risk_measure, cut_selection, parallel, regularisation, output)
    # catch InterruptException
        # warn("Terminating early")
    # end
    solution
end

terminate(converged::Bool, do_test::Bool) = converged & do_test
terminate(converged::Bool, do_test::Bool, simulation_passes::Range) = converged
terminate(converged::Bool, do_test::Bool, simulation_passes) = terminate(converged, do_test)

convergence_pass!(ty::Serial, m::SDDPModel, convergence, cut_selection::CutSelectionMethod) = convergence_pass!(m, convergence)
convergence_pass!(ty::ConvergenceTest, m::SDDPModel, convergence, cut_selection::CutSelectionMethod) = parallel_convergence_pass!(m, convergence, cut_selection)

backward_pass!(ty::BackwardPass, m::SDDPModel, method::CutSelectionMethod, regularisation::Regularisation) = parallel_backward_pass!(m, ty.cuts_per_processor, method)
backward_pass!(ty::Serial, m::SDDPModel, method::CutSelectionMethod, regularisation::Regularisation) = backward_pass!(m, method.frequency > 0, method, regularisation)

function solve!(m::SDDPModel, solution::Solution, convergence::Convergence, maximum_iterations::Int, risk_measure::RiskMeasure, cut_selection::CutSelectionMethod, parallel::Parallel, regularisation::Regularisation, output::Union{Void, ASCIIString})

    # Intialise model on worker processors
    if !isa(parallel, Parallel{Serial, Serial})
        nworkers = initialise_workers!(m)
    end

    # Initialise cut selection storage since we need it to communicate between processors
    initialise_cut_selection!(m)

    # Set risk aversion parameters
    m.beta_quantile, m.risk_lambda = risk_measure.beta, risk_measure.lambda

    iterations=0
    nsimulations = 0
    last_convergence_test, last_cutselection = 0, 0
    time_backwards, time_forwards, time_cutselection = 0., 0., 0.
    npass = isa(parallel.backward_pass, BackwardPass)?(parallel.backward_pass.cuts_per_processor * nworkers):1

    while iterations < maximum_iterations

        # Update iteration counters
        iterations += npass
        last_convergence_test += npass
        last_cutselection += npass

        # Cutting passes
        tic()
        backward_pass!(parallel.backward_pass, m, cut_selection, regularisation)
        time_backwards += toq()

        # Rebuild models if using Cut Selection
        if cut_selection.frequency > 0 && last_cutselection >= cut_selection.frequency
            tic()
            rebuild_stageproblems!(cut_selection, m)
            time_cutselection += toq()
        end

        # Test convergence if appropriate
        if convergence.frequency > 0 && last_convergence_test >= convergence.frequency
            last_convergence_test = 0
            tic()
            (is_converged, n) = convergence_pass!(parallel.convergence_test, m, convergence, cut_selection)
            nsimulations += Int(n)
            time_forwards += toq()

            # Output to user
            log!(solution, m.confidence_interval[1], m.confidence_interval[2], m.valid_bound, iterations, time_backwards, nsimulations, time_forwards, time_cutselection)
            print_stats(m, iterations, time_backwards, nsimulations, time_forwards, time_cutselection)
            if output != nothing
                open(output, "a") do outputfile
                    write(outputfile, textify(solution.trace[end]))
                end
            end

            # Terminate if converged and appropriate
            if terminate(is_converged, convergence.terminate, convergence.n)
                # We have terminated due to convergence test
                setStatus!(solution, CONVERGENCE_TERMINATION)
                return
            end
        end
    end
    # We have terminated due to iteration limit
    setStatus!(solution, ITERATION_TERMINATION)
end

function simulate(m::SDDPModel, n::Int, vars::Vector{Symbol}=[]; parallel=false)
    if parallel && length(workers()) < 2
        warn("Paralleisation requested but Julia is only running with a single processor. Start julia with `julia -p N` or use the `addprocs(N)` function to load N workers. Running in serial mode.")
        parallel = false
    end
    if parallel
        results = parallel_simulate(m, n, vars)
    else
        results = serial_simulate(m, n, vars)
    end
    results
end

function simulate(m::SDDPModel, vars::Vector{Symbol}=Symbol[]; markov=ones(Int, m.stages), kwargs...)
    n=1
    results = Dict{Symbol, Any}(:Objective=>zeros(Float64,n))
    for (s, t) in vcat(collect(zip(vars, fill(Any, length(vars)))), [(:Scenario, Int), (:Markov, Int), (:Future, Float64), (:Current, Float64)])
        results[s] = Array(Vector{t}, m.stages)
        for i=1:m.stages
            results[s][i] = Array(t, n)
        end
    end

    scenario = 0

    # Initialise the objective storage
    obj = 0.

    # For all stages
    for stage=1:m.stages
        # realise on scenario
        sp = m.stage_problems[stage, markov[stage]]
        scenario = load_scenario!(m, sp)
        data = stagedata(sp)
        for (key, series) in kwargs
            @assert haskey(data.scenario_constraint_names, key)
            cidx = data.scenario_constraint_names[key]
            JuMP.setRHS(data.scenario_constraints[cidx][1], series[stage])
        end

        # solve
        solve!(sp)

        # Add objective (stage profit only)
        obj += get_true_value(sp)

        # Save results if necesary
        store_results!(results, vars, sp, stage, 1, markov[stage], scenario)
        # pass forward if necessary
        if stage < m.stages
            pass_states_forward!(m, stage, markov[stage])
        end
    end

    results
end

end
