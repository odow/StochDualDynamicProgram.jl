#  Copyright 2016, Oscar Dowson

# Solve the classic hydro scheduling problem under price uncertainty.
#
# There are two reservoirs (upper and lower). In each time period, a reservoir can
#     release water through its turbine, spill water over the edge of the dam
#     (and therefore not through the turbine) and purchase water to top up the
#     s torage level.
#
# Each turbine has an identical piecewise linear response function.
#     f(water units) = electricity units

# Add processors
addprocs(4 - nprocs())

# We need to load out data everywhere
@everywhere begin
    using StochDualDynamicProgram, JuMP

    # For repeatability (note: don't use the same seed on two processors!)
    srand(11111 * myid())

    # Names of the reservoirs
    RESERVOIRS = [:upper, :lower]

    # Knots of the turbine response function
    A = [(50, 55), (60, 65), (70, 70)]

    # Prices[stage, markov state]
    Price = [
        1 2;
        2 1;
        3 4
    ]

    # Transition matrix
    Transition = [
        0.6 0.4;
        0.3 0.7
    ]

    # Price of purchasing and spilling water
    #   ($/Unit)
    M = 1000

    # Maximum level
    reservoir_max = Dict{Symbol, Float64}(
        :upper => 200,
        :lower => 200
        )

    # Initial fill
    reservoir_initial = Dict{Symbol, Float64}(
        :upper => 200,
        :lower => 200
        )

    # rainfal inflow
    rainfall = [
        10 20 30;
        20 30 40;
         0  5 10;
         0  5 10;
        10 15 20
    ]

    n = length(A)
end

# Initialise SDDP Model
m = SDDPModel(
        stages            = 3,
        scenarios         = 5,
        markov_states     = 2,
        transition        = Transition,
        value_to_go_bound = 1500
    ) do sp, stage, markov_state
    # ------------------------------------------------------------------
    #   SDDP State Variables
    # Level of upper reservoir
    @state(sp, 0 <= reservoir[r=RESERVOIRS] <= reservoir_max[r], reservoir0=reservoir_initial[r])
    # ------------------------------------------------------------------
    #   Additional variables
    @variables(sp, begin
        # Rainfall inflow into upper reservoir
        inflow >= 0

        # Quantity to flow through turbine of reservoir r
        outflow[r=RESERVOIRS] >= 0

        # Quantity to spill over edge of reservoir r
        spill[r=RESERVOIRS] >= 0

        # Total quantity of water
        generation_quantity >= 0

        # Proportion of levels to dispatch on
        0 <= dispatch[reservoir=RESERVOIRS, level=1:n] <= 1
    end)

    @constraints(sp, begin
        # ------------------------------------------------------------------
        # Conservation constraints
        reservoir[:upper] == reservoir0[:upper] - (outflow[:upper] + spill[:upper]) + inflow

        reservoir[:lower] == reservoir0[:lower] +
            (outflow[:upper] + spill[:upper]) -
            (outflow[:lower] + spill[:lower])

        # Dispatch combination of levels
        dispatched[reservoir=RESERVOIRS], sum{dispatch[reservoir, level], level=1:n} <= 1
    end)

    @relaxedconstraints(sp, 1000, begin
        # ------------------------------------------------------------------
        # Reservoir constraints
        # Flow out
        flowout[reservoir=RESERVOIRS], outflow[reservoir] == sum{A[level][1] * dispatch[reservoir, level], level=1:n}

        # Total quantity generated
        generation_quantity == sum{A[level][2] * dispatch[reservoir,level], reservoir=RESERVOIRS, level=1:n}
    end)

    # Random rainfall
    @scenarioconstraint(sp, scenario=1:5, inflow <= rainfall[scenario, stage])

    # ------------------------------------------------------------------
    #   Objective Function
    @stageprofit(sp, Price[stage, markov_state]*generation_quantity)
end

@time solvestatus = solve(m,
    maximum_iterations      = 50,
    expectedvalueiterations = 10,
    policy_estimation       = MonteCarloEstimator(
                                 frequency = 1,
                                 min       = 100,
                                 max       = 1000,
                                 step      = 100
                             ),
    forward_pass            = ForwardPass(
                                 scenarios = 10
                             ),
    backward_pass           = BackwardPass(
                                 multicut = true
                             ),
    parallel                = Parallel(),
    risk_measure            = NestedCVar(
                                 beta   = 0.5,
                                 lambda = 0.75
                             ),
    cut_selection           = LevelOne(5),
    print_level             = 2
)
@assert status(solvestatus) == :MaximumIterations

@time results = simulate(m,
    1000,
    [:reservoir, :dispatch, :outflow, :spill],
    parallel = true      # specify parallel simulate
    )

println("Mean of simulation was $(mean(results[:Objective]))")
