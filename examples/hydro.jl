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
    m = SDDPModel(stages=3, markov_states=2, scenarios=1, transition=Transition, value_to_go_bound=1500) do sp, stage, markov_state
        # ------------------------------------------------------------------
        #   SDDP State Variables
        # Level of upper reservoir
        @state(sp, 0 <= reservoir[r=RESERVOIRS] <= reservoir_max[r], reservoir0=reservoir_initial[r])
        # ------------------------------------------------------------------
        #   Additional variables
        @variables(sp, begin
            # Quantity to flow through turbine of reservoir r
            outflow[r=RESERVOIRS] >= 0

            # Quantity to spill over edge of reservoir r
            spill[r=RESERVOIRS] >= 0

            # Total quantity of water
            generation_quantity >= 0

            # Proportion of levels to dispatch on
            0 <= dispatch[reservoir=RESERVOIRS, level=1:n] <= 1
        end)

        # ------------------------------------------------------------------
        # Conservation constraints
        @constraints(sp, begin
            reservoir[:upper] == reservoir0[:upper] - (outflow[:upper] + spill[:upper])

            reservoir[:lower] == reservoir0[:lower] +
                (outflow[:upper] + spill[:upper]) -
                (outflow[:lower] + spill[:lower])
        end)

        # ------------------------------------------------------------------
        # Reservoir constraints
        for reservoir in RESERVOIRS
            @constraints(sp, begin
                # Flow out
                outflow[reservoir] == sum{
                    A[level][1] * dispatch[reservoir, level],
                    level=1:n}

                # Dispatch combination of levels
                sum{dispatch[reservoir, level], level=1:n} <= 1
            end)
        end

        # Total quantity generated
        @constraint(sp, generation_quantity == sum{
            A[level][2] * dispatch[reservoir,level],
            reservoir=RESERVOIRS, level=1:n}
        )

        # ------------------------------------------------------------------
        #   Objective Function
        @stageprofit(sp, Price[stage, markov_state]*generation_quantity)
    end

    m2 = copy(m)

    info("Sanity check. Perform more iterations than needed to check bounds do not cross.")
    @time solve(m,                # Solve the model using the SDDP algorithm
        convergence=Convergence(1000, 10),
        maximum_iterations=200
    )

    info("Cut selection comparison.")
    @time solve(m2,                # Solve the model using the SDDP algorithm
        convergence=Convergence(1000, 10),
        maximum_iterations=200,
        cut_selection = LevelOne(20)
    )

    results = simulate(m,   # Simulate the policy
        1000,               # number of monte carlo realisations
        [:reservoir,  # variables to return
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
        @state(sp,
            0 <= upper_reservoir <= reservoir_max[:upper],
            upper_reservoir0=reservoir_initial[:upper]
        )

        # Level of lower reservoir
        @state(sp,
            0 <= lower_reservoir <= reservoir_max[:lower],
            lower_reservoir0=reservoir_initial[:lower]
        )

        # ------------------------------------------------------------------
        #   Additional variables
        # Quantity to flow through turbine of reservoir r
        @variable(sp, outflow[r=RESERVOIRS] >= 0)

        # Quantity to spill over edge of reservoir r
        @variable(sp, spill[r=RESERVOIRS] >= 0)

        # Total quantity of water
        @variable(sp, generation_quantity >= 0)

        # Proportion of levels to dispatch on
        @variable(sp, 0 <= dispatch[reservoir=RESERVOIRS, level=1:n] <= 1)

        # ------------------------------------------------------------------
        #   Objective Function
        @stageprofit(sp, -Price[markov_state, stage]*generation_quantity)

        # ------------------------------------------------------------------
        # Conservation constraints
        @constraint(sp, upper_reservoir == upper_reservoir0 -
            (outflow[:upper] + spill[:upper])
        )

        @constraint(sp, lower_reservoir == lower_reservoir0 +
            (outflow[:upper] + spill[:upper]) -
            (outflow[:lower] + spill[:lower])
        )

        # ------------------------------------------------------------------
        # Reservoir constraints
        for reservoir in RESERVOIRS
            # Flow out
            @constraint(sp, outflow[reservoir] == sum{
                A[level][1] * dispatch[reservoir, level],
                level=1:n}
            )

            # Dispatch combination of levels
            @constraint(sp, sum{dispatch[reservoir, level], level=1:n} <= 1)
        end

        # Total quantity generated
        @constraint(sp, generation_quantity == sum{
            A[level][2] * dispatch[reservoir,level],
            reservoir=RESERVOIRS, level=1:n}
        )

    end

    solve(m,                # Solve the model using the SDDP algorithm
        convergence=Convergence(1000, 10),
        maximum_iterations=50
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
