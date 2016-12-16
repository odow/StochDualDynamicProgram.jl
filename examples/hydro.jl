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
using StochDualDynamicProgram, JuMP

# For repeatability
srand(11111)

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
        0 <= dispatch[r=RESERVOIRS, level=1:n] <= 1
    end)

    @constraints(sp, begin
        # ------------------------------------------------------------------
        # Conservation constraints
        reservoir[:upper] == reservoir0[:upper] - (outflow[:upper] + spill[:upper])

        reservoir[:lower] == reservoir0[:lower] +
            (outflow[:upper] + spill[:upper]) -
            (outflow[:lower] + spill[:lower])

        # Total quantity generated
        generation_quantity == sum(A[level][2] * dispatch[reservoir,level] for reservoir in RESERVOIRS for level in 1:n)

        # ------------------------------------------------------------------
        # Reservoir constraints
        # Flow out
        flowout[r=RESERVOIRS], outflow[r] == sum(A[level][1] * dispatch[r, level] for level in 1:n)

        # Dispatch combination of levels
        dispatched[r=RESERVOIRS], sum(dispatch[r, level] for level in 1:n) <= 1
    end)

    # ------------------------------------------------------------------
    #   Objective Function
    @stageprofit(sp, Price[stage, markov_state]*generation_quantity)
end

m2 = copy(m)

mcestimator = MonteCarloEstimator(
    frequency = 1,
    min       = 5,
    max       = 100,
    step      = 10
)

@time solvestatus = solve(m,
    maximum_iterations = 50,
    policy_estimation  = mcestimator,
    bound_convergence  = BoundConvergence(
                            after = 5,
                            tol   = 1e-10
                        ),
    print_level        = 2
)
@assert status(solvestatus) == :BoundConvergence

@time solvestatus = solve(m2,
    maximum_iterations = 20,
    policy_estimation  = mcestimator,
    forward_pass       = ForwardPass(
                            scenarios       = 1:10,
                            uniformsampling = true
                        ),
    cut_selection      = LevelOne(5),
    print_level        = 1
)
@assert status(solvestatus) == :MaximumIterations

results = simulate(m,   # Simulate the policy
    1000,               # number of monte carlo realisations
    [:reservoir,        # variables to return
    :dispatch,
    :outflow,
    :spill]
    )
