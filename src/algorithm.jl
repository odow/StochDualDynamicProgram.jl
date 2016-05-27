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
        # resize storage for forward pass
        resizeforwardstorage!(m, forward_scenarios)

        # forward pass
        forwardpass!(m, forward_scenarios, isa(cut_selection, LevelOne))

        # estimate bound
        setCI!(m, estimatebound(getobj(m), 0.95))
        log.ci_lower, log.ci_upper = m.confidence_interval

        # backward pass
        backwardpass!(m, risk_measure, regularisation)
        log.cuts += forward_scenarios

        # Calculate a new upper bound
        setbound!(m)
        log.bound = m.valid_bound

        # run cut selection
        cutselection!(m, cut_selection, solution.iterations)

        solution.iterations += forward_scenarios
        push!(solution.trace, copy(log))
        print(m, log)

        if cut_output_file != nothing
            writecuts!(m, cut_output_file, cutswrittentofile)
            cutswrittentofile += forward_scenarios
        end
    end
end

function cutselection!(m::SDDPModel, cutselection::CutSelectionMethod, iteration)
    if mod(iteration, cutselection.frequency) == 0
        println("Running cut selection")
        rebuild_stageproblems!(m, cutselection)
    end
end
cutselection!(m::SDDPModel, cutselection::NoSelection, iteration) = nothing

function setriskmeasure!(m::SDDPModel, riskmeasure::RiskMeasure)
    for sp in m.stage_problems
        stagedata(sp).beta_quantile = riskmeasure.beta
        stagedata(sp).lambda_weight = riskmeasure.lambda
    end
end

getobj(m::SDDPModel) = m.forwardstorage.obj
function estimatebound(obj::Vector{Float64}, conflevel)
    if length(obj) > 5
        return confidenceinterval(obj, conflevel)
    else
        return (mean(obj), mean(obj))
    end
end

function confidenceinterval(x, conf_level=0.95)
    tstar = quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    SE = std(x)/sqrt(length(x))
    lo, hi = mean(x) + [-1, 1] * tstar * SE
    return (lo, hi)
end
setCI!{T<:Real}(m::SDDPModel, v::Tuple{T,T}) = (m.confidence_interval = v)

function resizeforwardstorage!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, n::Int)
    resize!(m.forwardstorage.obj, n)
    nx = length(stagedata(m, 1,1).state_vars)
    for i=(m.forwardstorage.n+1):n
        push!(m.forwardstorage.x, zeros(nx, T))
        push!(m.forwardstorage.W, zeros(Int, T))
    end
end

function getstagevalue(sp::Model)
    @assert issubproblem(sp)
    if stagedata(sp).theta != nothing
        return (getobjectivevalue(sp) - getvalue(stagedata(sp).theta))::Float64
    else
        return getobjectivevalue(sp)::Float64
    end
end

"""
Calculate the upper bound of the first stage problem
"""
function setbound!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM})
    # Initialise
    obj = 0.0
    if m.initial_markov_state == 0
        # Lets average over all first stage probles (markov states x scenarios)
        for i=1:M
            sp = subproblem(m, 1, i)        # get subproblem
            for s=1:S
                load_scenario!(sp, s)   # realise scenario
                forwardsolve!(sp)          # solve
                obj += transitionprobability(m, 1, m.initial_markov_state, i) * m.scenario_probability[s] * getobjectivevalue(sp)
            end
        end
    else
        # Lets just  average over the scenarios in the initial markov state
        sp = subproblem(m, 1, m.initial_markov_state) # get subproblem
        for s=1:S
            load_scenario!(sp, s)                  # load scenario
            forwardsolve!(sp)                         # solve
            obj += m.scenario_probability[s] * getobjectivevalue(sp)
        end
    end
    setbound!(m, obj)    # Update the bound
end
