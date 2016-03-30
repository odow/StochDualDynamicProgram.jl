using StochDualDynamicProgram, JuMP

# For repeatability
srand(11111)


# Names of the reservoirs
RESERVOIRS = [:upper, :lower]

# Water duty converts cu m to MWh
W = Dict{Symbol, Float64}(
:upper => 1.0,
:lower => 1.0
)

# Prices[stage]
Price = [
    1 ;
    2 ;
    3 ;
	4 ;
	5
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

# Initialise SDDP Model
m = SDDPModel(stages=5, value_to_go_bound=10000) do sp, stage
    # ------------------------------------------------------------------
    #   SDDP State Variables
    # Level of upper reservoir
    @defStateVar(sp, 0 <= reservoir[r=RESERVOIRS] <= reservoir_max[r], reservoir0=reservoir_initial[r])

    # ------------------------------------------------------------------
    #   Additional variables
    # Quantity to flow through turbine of reservoir r
    @defVar(sp, outflow[r=RESERVOIRS] >= 0)

    # Quantity to spill over edge of reservoir r
    @defVar(sp, spill[r=RESERVOIRS] >= 0)

	# Generation of reservoir r
    @defVar(sp, generation[r=RESERVOIRS] >= 0)

    # Total quantity of water
    @defVar(sp, generation_quantity >= 0)


    # ------------------------------------------------------------------
    # Conservation constraints
    @addConstraint(sp, reservoir[:upper] == reservoir0[:upper] -
        (outflow[:upper] + spill[:upper])
    )

    @addConstraint(sp, reservoir[:lower] == reservoir0[:lower] +
        (outflow[:upper] + spill[:upper]) -
        (outflow[:lower] + spill[:lower])
    )

    # ------------------------------------------------------------------
    # Reservoir constraints
    for r in RESERVOIRS
        # Flow out
        @addConstraint(sp, generation[r] == W[r] * outflow[r]
        )

    end

    # Total quantity generated
    @addConstraint(sp, generation_quantity == sum{
        generation[r],
        r=RESERVOIRS}
    )

    # ------------------------------------------------------------------
    #   Objective Function
    @setStageProfit(sp, Price[stage]*generation_quantity)

end

solve(m,                # Solve the model using the SDDP algorithm
    simulation_passes=1000,
    log_frequency=1,
    maximum_iterations=5
)

results = simulate(m,   # Simulate the policy
    1,               # number of monte carlo realisations
	[:reservoir,  # variables to return
    :generation,
    :outflow,
    :spill]
    )
println("Optimal Solution is:")
for stage=1:5
    println("At stage ", stage)
    println(results[:outflow][stage][1])
#	println(results[:generation][stage][1])
    println(results[:reservoir][stage][1])
    println(" ")
end

results[:Objective][1]
