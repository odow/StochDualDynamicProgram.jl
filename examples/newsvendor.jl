using StochDualDynamicProgram, JuMP

function solve_newsvendor()
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
            stages=3,
            markov_states=2,
            scenarios=2,
            transition=Transition,
            initial_markov_state=1,
            value_to_go_bound = 1000,
            cuts_filename="news_vendor.csv"
        ) do sp, stage, markov_state

        # ====================
        #   State variable
        @defStateVar(sp, 0 <= stock <= 100, stock0=5)

        # ====================
        #   Other variables
        @defVar(sp, buy  >= 0)  # Quantity to buy
        @defVar(sp, sell >= 0)  # Quantity to sell

        # ====================
        #   Scenarios
        @addScenarioConstraint(sp, D=Demand[stage,:], sell <= D)

        # ====================
        #   Objective
        @setStageProfit(sp, sell * RetailPrice - buy * PurchasePrice[markov_state])

        # ====================
        #   Dynamics constraint
        @addConstraint(sp, stock == stock0 + buy - sell)

    end

    # load_cuts!(m)

    solve(m,                # Solve the model using the SDDP algorithm
        simulation_passes=1000,
        log_frequency=10,
        maximum_iterations=50,
        beta_quantile=0.6,
        risk_lambda = 0.5
    )

    results = simulate(m,   # Simulate the policy
        1000,               # number of monte carlo realisations
        [:stock, :buy, :sell]
        )

    return results
end

function solve_newsvendor2()
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
            stages=3,
            markov_states=2,
            scenarios=2,
            transition=Transition,
            initial_markov_state=1,
            value_to_go_bound = 1000,
        ) do sp, stage, markov_state

        # ====================
        #   State variable
        @defStateVar(sp, 0 <= stock <= 100, stock0=5)

        # ====================
        #   Other variables
        @defVar(sp, buy  >= 0)  # Quantity to buy
        @defVar(sp, sell >= 0)  # Quantity to sell

        # ====================
        #   Scenarios
        @addScenarioConstraint(sp, D=Demand[stage,:], sell <= D)

        # ====================
        #   Objective
        @setStageProfit(sp, sell * RetailPrice - buy * PurchasePrice[markov_state])

        # ====================
        #   Dynamics constraint
        @addConstraint(sp, stock == stock0 + buy - sell)

    end

    m1 = copy(m)
    m2 = copy(m)

    srand(11111)
    info("Don't check for duplicate cuts")
    @time solve(m,                # Solve the model using the SDDP algorithm
        simulation_passes=10000,
        log_frequency=30,
        maximum_iterations=30
    )

    srand(11111)
    info("Check for duplicate cuts")
    @time solve(m1,                # Solve the model using the SDDP algorithm
        simulation_passes=10000,
        log_frequency=30,
        maximum_iterations=30,
        check_duplicate_cuts=true
    )

    srand(11111)
    info("Solving using varying number of simulation passes")
    @time solve(m2,                # Solve the model using the SDDP algorithm
        simulation_passes=linspace(100, 10000, 10),
        log_frequency=1,
        maximum_iterations=50
    )

    results = simulate(m,   # Simulate the policy
        1000,               # number of monte carlo realisations
        [:stock, :buy, :sell]
        )

    return results
end
