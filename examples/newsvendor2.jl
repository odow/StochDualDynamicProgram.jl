#  Copyright 2016, Oscar Dowson

using StochDualDynamicProgram, JuMP

srand(11111)

# Demand for newspapers
# There are two equally probable scenarios in each stage
#   Demand[stage, scenario]
Demand = [
    10. 15.;
    12. 20.;
    8.  20.
]

# Markov state purchase prices
PurchasePrice = [5., 8.]

RetailPrice = 7.

# Transition matrix
Transition = Array{Float64, 2}[
    [0.6 0.4; 0.3 0.7],
    [0.3 0.7; 0.3 0.7],
    [0.5 0.5; 0.5 0.5]
  ]

# Initialise SDDP Model
m = SDDPModel(
        stages               = 3,
        markov_states        = 2,
        scenarios            = 2,
        transition           = Transition,
        initial_markov_state = 1,
        value_to_go_bound    = 1000
    ) do sp, stage, markov_state

    # ====================
    #   State variable
    @state(sp, 0 <= stock <= 100, stock0=5)

    # ====================
    #   Other variables
    @variables(sp, begin
        buy  >= 0  # Quantity to buy
        sell >= 0  # Quantity to sell
    end)

    # ====================
    #   Scenarios
    @scenarioconstraint(sp, scenario=1:2, sell <= Demand[stage,scenario])

    # ====================
    #   Objective
    @stageprofit(sp, sell * RetailPrice - buy * PurchasePrice[markov_state])

    # ====================
    #   Dynamics constraint
    @constraint(sp, stock == stock0 + buy - sell)

end

m1 = copy(m)
m2 = copy(m)

montecarlo = MonteCarloEstimator(
    frequency = 10,
    min       = 100,
    max       = 1000,
    step      = 100,
    terminate = true
)

@time solvestatus = solve(m,
    maximum_iterations      = 30,
    expectedvalueiterations = 10,
    policy_estimation       = montecarlo,
    cut_output_file         = "news_vendor.csv"
)
@assert status(solvestatus) == :PolicyConverence

@time solvestatus = solve(m1,
    maximum_iterations = 30,
    bound_convergence  = BoundConvergence(
                            after = 5,
                            tol   = 1e-5
                        ),
    cut_selection      = LevelOne(5)
)
@assert status(solvestatus) == :BoundConvergence

results = simulate(m,   # Simulate the policy
    1000,               # number of monte carlo realisations
    [:stock, :buy, :sell]
    )

info("Rebuilding model from cuts")

# Initialise SDDP Model
m2 = SDDPModel(
        stages=3,
        markov_states=2,
        scenarios=2,
        transition=Transition,
        initial_markov_state=1,
        value_to_go_bound = 1000,
    ) do sp, stage, markov_state

    @state(sp, 0 <= stock <= 100, stock0=5)

    @variables(sp, begin
        buy  >= 0  # Quantity to buy
        sell >= 0  # Quantity to sell
    end)

    @scenarioconstraint(sp, scenario=1:2, sell <= Demand[stage,scenario])

    @stageprofit(sp, sell * RetailPrice - buy * PurchasePrice[markov_state])

    @constraint(sp, stock == stock0 + buy - sell)
end

loadcuts!(m2, "news_vendor.csv")

results1 = simulate(m2,   # Simulate the policy
    1000,               # number of monte carlo realisations
    [:stock, :buy, :sell]
    )

# Hopefully the mean of the simulated in memory is within $1 of the simulated
# from cut files
@assert abs(mean(results[:Objective]) - mean(results1[:Objective])) < 1

rm("news_vendor.csv")
