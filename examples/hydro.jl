using StochDualDynamicProgram, JuMP

# For repeatability
srand(11111)

"""
Solve the classic hydro scheduling problem under price uncertainty.

There are two reservoirs (upper and lower). In each time period, a reservoir can
    release water through its turbine, spill water over the edge of the dam
    (and therefore not through the turbine) and purchase water to top up the
    storage level.

Each turbine has an identical piecewise linear response function.
    f(water units) = electricity units
"""
function solve_hydro()
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

    n = length(A)

    # Initialise SDDP Model
    m = SDDPModel(stages=3, markov_states=2, scenarios=1, transition=Transition) do sp, stage, markov_state
        # ------------------------------------------------------------------
        #   SDDP State Variables
        # Level of upper reservoir
        @defStateVar(sp,
            0 <= upper_reservoir <= reservoir_max[:upper],
            upper_reservoir0==reservoir_initial[:upper]
        )

        # Level of lower reservoir
        @defStateVar(sp,
            0 <= lower_reservoir <= reservoir_max[:lower],
            lower_reservoir0==reservoir_initial[:lower]
        )

        # ------------------------------------------------------------------
        #   SDDP Value to Go
        @defValueToGo(sp, value_to_go <= M)
        if stage==3
            # No value to go in last stage
            @addConstraint(sp, value_to_go==0)
        end
        # ------------------------------------------------------------------
        #   Additional variables
        # Quantity to flow through turbine of reservoir r
        @defVar(sp, outflow[r=RESERVOIRS] >= 0)

        # Quantity to spill over edge of reservoir r
        @defVar(sp, spill[r=RESERVOIRS] >= 0)

        # Total quantity of water
        @defVar(sp, generation_quantity >= 0)

        # Proportion of levels to dispatch on
        @defVar(sp, 0 <= dispatch[reservoir=RESERVOIRS, level=1:n] <= 1)

        # ------------------------------------------------------------------
        # Conservation constraints
        @addConstraint(sp, upper_reservoir == upper_reservoir0 -
            (outflow[:upper] + spill[:upper])
        )

        @addConstraint(sp, lower_reservoir == lower_reservoir0 +
            (outflow[:upper] + spill[:upper]) -
            (outflow[:lower] + spill[:lower])
        )

        # ------------------------------------------------------------------
        # Reservoir constraints
        for reservoir in RESERVOIRS
            # Flow out
            @addConstraint(sp, outflow[reservoir] == sum{
                A[level][1] * dispatch[reservoir, level],
                level=1:n}
            )

            # Dispatch combination of levels
            @addConstraint(sp, sum{dispatch[reservoir, level], level=1:n} <= 1)
        end

        # Total quantity generated
        @addConstraint(sp, generation_quantity == sum{
            A[level][2] * dispatch[reservoir,level],
            reservoir=RESERVOIRS, level=1:n}
        )

        # ------------------------------------------------------------------
        #   Objective Function
        @setObjective(sp, Max,
            Price[stage, markov_state]*generation_quantity +
            value_to_go
        )
    end

    solve(m,                # Solve the model using the SDDP algorithm
        forward_passes=1000,  # number of realisations in bound simulation
        backward_passes=10,   # number of cutting iterations before convergence check
        max_iterations=5
    )

    results = simulate(m,   # Simulate the policy
        1000,               # number of monte carlo realisations
        [:lower_reservoir,  # variables to return
        :upper_reservoir,
        :dispatch,
        :outflow,
        :spill]
        )

    return results
end

# Hydro problem with different transition matrix
function solve_hydro2()
    # Names of the reservoirs
    RESERVOIRS = [:upper, :lower]

    # Knots of the turbine response function
    A = [(50, 55), (60, 65), (70, 70)]


    # Prices[markov state, scenario]
    Price = [
        5 5.5 6;    # high price
        1 2 3       # low price
    ]

    # Transition matrix
    Transition = Array{Float64, 2}[
        [0.6 0.4;0.3 0.7],
        [0.3 0.7;0.3 0.7],
        [0.5 0.5;0.5 0.5]
      ]

    # Price of purchasing and spilling water
    #   ($/Unit)
    M = 2000

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

    n = length(A)

    # Initialise SDDP Model
    m = SDDPModel(sense=:Min, stages=3, markov_states=2, transition=Transition, value_to_go_bound=-1500) do sp, stage, markov_state
        # ------------------------------------------------------------------
        #   SDDP State Variables
        # Level of upper reservoir
        @defStateVar(sp,
            0 <= upper_reservoir <= reservoir_max[:upper],
            upper_reservoir0==reservoir_initial[:upper]
        )

        # Level of lower reservoir
        @defStateVar(sp,
            0 <= lower_reservoir <= reservoir_max[:lower],
            lower_reservoir0==reservoir_initial[:lower]
        )

        # ------------------------------------------------------------------
        #   Additional variables
        # Quantity to flow through turbine of reservoir r
        @defVar(sp, outflow[r=RESERVOIRS] >= 0)

        # Quantity to spill over edge of reservoir r
        @defVar(sp, spill[r=RESERVOIRS] >= 0)

        # Total quantity of water
        @defVar(sp, generation_quantity >= 0)

        # Proportion of levels to dispatch on
        @defVar(sp, 0 <= dispatch[reservoir=RESERVOIRS, level=1:n] <= 1)

        # ------------------------------------------------------------------
        #   Objective Function
        @setStageProfit(sp, -Price[markov_state, stage]*generation_quantity)

        # ------------------------------------------------------------------
        # Conservation constraints
        @addConstraint(sp, upper_reservoir == upper_reservoir0 -
            (outflow[:upper] + spill[:upper])
        )

        @addConstraint(sp, lower_reservoir == lower_reservoir0 +
            (outflow[:upper] + spill[:upper]) -
            (outflow[:lower] + spill[:lower])
        )

        # ------------------------------------------------------------------
        # Reservoir constraints
        for reservoir in RESERVOIRS
            # Flow out
            @addConstraint(sp, outflow[reservoir] == sum{
                A[level][1] * dispatch[reservoir, level],
                level=1:n}
            )

            # Dispatch combination of levels
            @addConstraint(sp, sum{dispatch[reservoir, level], level=1:n} <= 1)
        end

        # Total quantity generated
        @addConstraint(sp, generation_quantity == sum{
            A[level][2] * dispatch[reservoir,level],
            reservoir=RESERVOIRS, level=1:n}
        )

    end

    solve(m,                # Solve the model using the SDDP algorithm
        forward_passes=2000,  # number of realisations in bound simulation
        backward_passes=10,   # number of cutting iterations before convergence check
        max_iterations=5
    )

    results = simulate(m,   # Simulate the policy
        1000,               # number of monte carlo realisations
        [:lower_reservoir,  # variables to return
        :upper_reservoir,
        :dispatch,
        :outflow,
        :spill]
        )

    return results
end
