# Solve the classic hydro scheduling problem under price uncertainty.
#
# There are two reservoirs (upper and lower). In each time period, a reservoir can
#     release water through its turbine, spill water over the edge of the dam
#     (and therefore not through the turbine) and purchase water to top up the
#     s torage level.
#
# Each turbine has an identical piecewise linear response function.
#     f(water units) = electricity units
using StochDualDynamicProgram, JuMP, CPLEX

# For repeatability
srand(11111)

# Names of the reservoirs
RESERVOIRS = [:upper, :lower]

# Knots of the turbine response function
A = [(0,0), (50, 55), (60, 65), (70, 70)]

# Prices[markov state, stage]
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
m = SDDPModel(
                sense = :Min,
               stages = 3,
        markov_states = 2,
           transition = Transition,
    value_to_go_bound = -1500,
               solver = CplexSolver(CPX_PARAM_SCRIND=0)
                                                   ) do sp, stage, markov_state
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

        power[r=RESERVOIRS]
    end)

    @constraints(sp, begin
        # ------------------------------------------------------------------
        # Conservation constraints
        reservoir[:upper] == reservoir0[:upper] - (outflow[:upper] + spill[:upper])

        reservoir[:lower] == reservoir0[:lower] +
            (outflow[:upper] + spill[:upper]) -
            (outflow[:lower] + spill[:lower])

        # Total quantity generated
        generation_quantity == sum{power[r], r=RESERVOIRS}
        powercurve[r=RESERVOIRS], power[r] == SOSII!(sp, outflow[r], A)
    end)

    # ------------------------------------------------------------------
    #   Objective Function
    @stageprofit(sp, -Price[markov_state, stage]*generation_quantity)
end

@time solvestatus = solve(m,
    maximum_iterations = 50,
    policy_estimation  = MonteCarloEstimator(
                            frequency          = 1,
                            min                = 5,
                            max                = 100,
                            step               = 10,
                            terminate          = false,
                            antitheticvariates = false
                        ),
    print_level        = 4
)
@assert status(solvestatus) == :MaximumIterations

results = simulate(m,   # Simulate the policy
    1000,               # number of monte carlo realisations
    [:reservoir]  # variables to return
    )
