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
    @scenarioconstraints(sp, D=Demand[stage,:], begin
        maxdemand, sell <= D
        mindemand, sell >= 0.5D
    end)

    # ====================
    #   Objective
    @stageprofit(sp, sell * RetailPrice - buy * PurchasePrice[markov_state])

    # ====================
    #   Dynamics constraint
    @constraint(sp, stock == stock0 + buy - sell)

end

@time solvestatus = solve(m,
    maximum_iterations = 50,
    policy_estimation  = MonteCarloEstimator(
                            frequency = 10,
                            min       = 5,
                            max       = 50,
                            step      = 5,
                            terminate = false
                        ),
    forward_pass       = ForwardPass(
                            regularisation = LinearRegularisation(1., 0.95)
                        ),
    risk_measure       = NestedCVar(
                            beta   = 0.6,
                            lambda = 0.5
                        )
)
@assert status(solvestatus) == :MaximumIterations

# Historical simulation
results = historicalsimulation(m,   # Simulate the policy
    [:stock, :buy, :sell],
    markov = [1, 2, 1],
    maxdemand = [10, 10, 10]
    )
@assert abs(results[:Objective][1] - 85) < 1e-5

# Monte-carlo simulation
results = simulate(m,   # Simulate the policy
    1000,               # number of monte carlo realisations
    [:stock, :buy, :sell]
    )
