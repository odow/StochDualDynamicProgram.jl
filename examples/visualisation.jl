# This example demonstates the D3.js visualisation
using JuMP, StochDualDynamicProgram

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

    @constraints(sp, begin
        # ------------------------------------------------------------------
        # Conservation constraints
        reservoir[:upper] == reservoir0[:upper] - (outflow[:upper] + spill[:upper])

        reservoir[:lower] == reservoir0[:lower] +
            (outflow[:upper] + spill[:upper]) -
            (outflow[:lower] + spill[:lower])

        # Total quantity generated
        generation_quantity == sum{A[level][2] * dispatch[reservoir,level], reservoir=RESERVOIRS, level=1:n}

        # ------------------------------------------------------------------
        # Reservoir constraints
        # Flow out
        flowout[reservoir=RESERVOIRS], outflow[reservoir] == sum{A[level][1] * dispatch[reservoir, level], level=1:n}

        # Dispatch combination of levels
        dispatched[reservoir=RESERVOIRS], sum{dispatch[reservoir, level], level=1:n} <= 1
    end)

    # ------------------------------------------------------------------
    #   Objective Function
    @stageprofit(sp, Price[stage, markov_state]*generation_quantity)
end

@time solve(m,                # Solve the model using the SDDP algorithm
    maximum_iterations=20
)

results = simulate(m,   # Simulate the policy
    100,               # number of monte carlo realisations
    [:reservoir]  # variables to return
    )

@visualise(results, (stage, replication), begin
	results[:Current][stage][replication],              (title="Accumulated Profit", ylabel="Accumulated Profit (\$)", cumulative=true)
	results[:Current][stage][replication],              (title="Weekly Income",      ylabel="Week Profit (\$)")
	results[:reservoir][stage][replication][:upper],    (title="Upper Reservoir",    ylabel="Level")
	results[:reservoir][stage][replication][:lower],    (title="Lower Reservoir")
	Price[stage, results[:Markov][stage][replication]], (ylabel="Price")
    results[:Future][stage][replication]
end)
