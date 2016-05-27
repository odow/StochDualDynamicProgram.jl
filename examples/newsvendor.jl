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
            value_to_go_bound = 1000
        ) do sp, stage, markov_state

        # ====================
        #   State variable
        @state(sp, 0 <= stock <= 100, stock0=5)

        # ====================
        #   Other variables
        @variable(sp, buy  >= 0)  # Quantity to buy
        @variable(sp, sell >= 0)  # Quantity to sell

        # ====================
        #   Scenarios
        @scenarioconstraint(sp, demand, D=Demand[stage,:], sell <= D)

        # ====================
        #   Objective
        @stageprofit(sp, sell * RetailPrice - buy * PurchasePrice[markov_state])

        # ====================
        #   Dynamics constraint
        @constraint(sp, stock == stock0 + buy - sell)

    end

    solve(m,                # Solve the model using the SDDP algorithm
        convergence=Convergence(1000, 10),
        maximum_iterations=50,
        risk_measure = NestedCVar(beta=0.6, lambda=0.5)
    )

    # # Historical simulation
    # results = simulate(m,   # Simulate the policy
    #     [:stock, :buy, :sell],
    #     demand=[10, 10, 10]
    #     )

    # Monte-carlo simulation
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
        @state(sp, 0 <= stock <= 100, stock0=5)

        # ====================
        #   Other variables
        @variable(sp, buy  >= 0)  # Quantity to buy
        @variable(sp, sell >= 0)  # Quantity to sell

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

    srand(11111)
    info("Don't check for duplicate cuts")
    @time solve(m,                # Solve the model using the SDDP algorithm
        convergence=Convergence(10000, 30),
        maximum_iterations=30,
        cut_output_file = "news_vendor.csv"
    )

    srand(11111)
    info("Cut selection")
    @time solve(m1,                # Solve the model using the SDDP algorithm
        convergence=Convergence(10000, 30),
        maximum_iterations=30,
        cut_selection = LevelOne(5)
    )

    srand(11111)
    info("Solving using varying number of simulation passes")
    @time solve(m2,                # Solve the model using the SDDP algorithm
        convergence=Convergence(linspace(100, 10000, 10), 1),
        maximum_iterations=50
    )

    results = simulate(m,   # Simulate the policy
        1000,               # number of monte carlo realisations
        [:stock, :buy, :sell]
        )

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

    @show mean(results[:Objective]), mean(results1[:Objective])

    rm("news_vendor.csv")

    return results
end
