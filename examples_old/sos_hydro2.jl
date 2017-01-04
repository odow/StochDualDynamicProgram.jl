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
using StochDualDynamicProgram, JuMP, CPLEX

# For repeatability
srand(11111)

# Names of the reservoirs, stations into a network
RESERVOIRS = Symbol[:Res1,:Res2]
STATIONS = Symbol[:Stn1,:Stn2]
SEA = Symbol[:Sea]
NODES = union(RESERVOIRS,STATIONS,SEA)

# HEADS dictionary maps the reservoirs with the stations for which
# head level is modeled
HEADS = Dict(
    :Stn1 => [:Res1],
    :Stn2 => [:Res2]
)

# Upstream dictionary maps the nodes Upstream
UPSTREAM_OUTFLOW = Dict(
    :Res1 => [],
    :Res2 => [:Stn1],
    :Stn1 => [:Res1],
    :Stn2 => [:Res2],
    :Sea => [:Stn2]
)

UPSTREAM_SPILL = Dict(
    :Res1 => [],
    :Res2 => [:Res1],
    :Stn1 => [],
    :Stn2 => [],
    :Sea => [:Res2]
)

TURBINES = Dict(
  :Stn1 => 2,
  :Stn2 => 3
)

# Knots of the turbine response function
TurbineCurves = Dict(
  :Stn1 => [(0,0), (50, 55), (60, 65), (70, 70)],
  :Stn2 => [(0,0), (50, 55), (60, 65), (70, 70)]
)

# PriceStates[markov state]
PriceStates = [1, 2, 3]

# Transition matrix
Transition = [
    0.6 0.4;
    0.3 0.7
]

# Price of purchasing and spilling water
#   ($/Unit)
M = 1000

# Maximum level
reservoir_min = Dict{Symbol, Any}(
    :Res1 => 0,
    :Res2 => 0
)

reservoir_max = Dict{Symbol, Any}(
    :Res1 => 1000,
    :Res2 => 1000
)

# Initial fill
reservoir_initial = Dict{Symbol, Any}(
    :Res1 => 200,
    :Res2 => 300
    )

#Flow bounds
spill_max = Dict(
    :Res1 =>70,
    :Res2 =>70
)

turbineFlow_min = Dict(
    :Stn1 => 0,
    :Stn2 => 0
)

turbineFlow_max = Dict(
    :Stn1 => 70,
    :Stn2 => 70
)

turbinePower_min = Dict(
    :Stn1 => 0,
    :Stn2 => 0
)

turbinePower_max = Dict(
    :Stn1 => 70,
    :Stn2 => 70
)

# Initialise SDDP Model
m = SDDPModel(
                sense = :Max,
               stages = 3,
        markov_states = 2,
            scenarios = 1,
           transition = Transition,
    value_to_go_bound = 100000,
               solver = CplexSolver(CPX_PARAM_SCRIND=0)
                                                    ) do sp, stage, markov_state
    # ------------------------------------------------------------------
    #   SDDP State Variables
    # Level of upper reservoir
    @state(sp, reservoir_min[r] <= reservoir[r=RESERVOIRS] <= reservoir_max[r], reservoir0=reservoir_initial[r])

    # ------------------------------------------------------------------
    #   Additional variables
    @variables(sp, begin
        # Quantity to flow through turbine of reservoir r
        outflow[n=NODES] >= 0

        # Quantity to spill over edge of reservoir r
        spill[r=RESERVOIRS] >= 0

        #Discharge quantity from the station
        stationFlow[s=STATIONS] >= 0
        turbineFlow[s=STATIONS,turbine=1:TURBINES[s]] >= 0

        #Power quantity from the station
        stationPower[s=STATIONS] >= 0
        turbinePower[s=STATIONS,turbine=1:TURBINES[s]] >= 0

        # Total quantity of water
        OfferQuantity >= 0

        # Number of units on
        0  <= unitsOnBool[s=STATIONS,turbine=1:TURBINES[s]]  <= 1, Bin
    end)

    @constraints(sp, begin

        # ------------------------------------------------------------------
        # Offer constraints
        OfferQuantity == sum{stationPower[s],s=STATIONS}

        # ------------------------------------------------------------------
        # Number of units on bool upper
        # numberUnitsOnUB[s=STATIONS], sum{unitsOnBool[s,turbine],turbine=1:TURBINES[s]} <= 1
        numberUnitsOnUB[s=STATIONS, t=2:TURBINES[s]], unitsOnBool[s,t] <= unitsOnBool[s,t-1]

        # ------------------------------------------------------------------
        # Station power definition
        StationPowerDefinition[s=STATIONS], stationPower[s] <= sum{turbinePower[s, turbine],turbine=1:TURBINES[s]}

        # ------------------------------------------------------------------
        # Station power definition
        stationFlowDefinition[s=STATIONS], stationFlow[s] == sum{turbineFlow[s, turbine],turbine=1:TURBINES[s]}

        # ------------------------------------------------------------------
        # Turbine Flow Bounds
        turbineFlowLB[s=STATIONS,t=1:TURBINES[s]], turbineFlow[s,t] >= turbineFlow_min[s]*unitsOnBool[s,t]
        turbineFlowUB[s=STATIONS,t=1:TURBINES[s]], turbineFlow[s,t] <= turbineFlow_max[s]*unitsOnBool[s,t]

        # ------------------------------------------------------------------
        # Turbine Flow SOSII
        turbinesos[s=STATIONS,t=1:TURBINES[s]], turbinePower[s, t] <= SOSII!(sp, turbineFlow[s,t], TurbineCurves[s])

        # ------------------------------------------------------------------
        # Water balance constraints
        WaterBalance[r=RESERVOIRS], reservoir[r] == reservoir0[r] + sum{outflow[n],n=UPSTREAM_OUTFLOW[r]} +
            sum{spill[n],n=UPSTREAM_SPILL[r]} - outflow[r] - spill[r]

        # ------------------------------------------------------------------
        # station flow balance constraints
        StationBalanceIn[s=STATIONS], stationFlow[s] == sum{outflow[n],n=UPSTREAM_OUTFLOW[s]} +
                                                        sum{spill[n],n=UPSTREAM_SPILL[s]}

        StationBalanceOUt[s=STATIONS], stationFlow[s] == outflow[s]

        # ------------------------------------------------------------------
        # Spill upper bound constraints
        SpillUB[r=RESERVOIRS], spill[r]<=spill_max[r]

    end)

    # ------------------------------------------------------------------
    #   Objective Function
    @stageprofit(sp, PriceStates[markov_state]*OfferQuantity)
end

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

results = simulate(m, 100, [:reservoir, :unitsOnBool, :stationPower, :stationFlow])

@visualise(results, (t,replication), begin
    results[:reservoir][t][replication][:Res1], (ylabel="Storage", title="Reservoir 1")
    results[:reservoir][t][replication][:Res2], (ylabel="Storage", title="Reservoir 2")
    sum([results[:unitsOnBool][t][replication][:Stn1, tur] for tur=1:TURBINES[:Stn1]]), (title="Unit Committment 1")
    sum([results[:unitsOnBool][t][replication][:Stn2, tur] for tur=1:TURBINES[:Stn2]]), (title="Unit Committment 2")
    results[:stationPower][t][replication][:Stn1], (title="Station Power 1")
    results[:stationPower][t][replication][:Stn2], (title="Station Power 2")
    results[:stationFlow][t][replication][:Stn1], (title="Station Flow 1")
    results[:stationFlow][t][replication][:Stn2], (title="Station Flow 2")
end)
