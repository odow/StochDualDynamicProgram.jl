# This makes a pretty picture showing how difference levels of risk aversion skew the profit distribution
using StochDualDynamicProgram, JuMP, Gadfly, DataFrames

# Hydro problem with different transition matrix
function solve_newsvendor(Demand, beta_quant=0.5, lambda=1.)
    # Demand for newspapers
    # There are two equally probable scenarios in each stage
    #   Demand[stage, scenario]
    # Demand = [
    #     10. 15.;
    #     12. 20.;
    #     8.  20.
    # ]

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
            scenarios=size(Demand)[2],
            transition=Transition,
            initial_markov_state=1,
            value_to_go_bound = 1000
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

    solve(m,                # Solve the model using the SDDP algorithm
        simulation_passes=1000,
        log_frequency=10,
        maximum_iterations=50,
        beta_quantile=beta_quant,
        risk_lambda = lambda
    )

    results = simulate(m,   # Simulate the policy
        3000,               # number of monte carlo realisations
        [:stock, :buy, :sell]
        )

    return results
end

Demand = 20 * rand(3,30)

r1 = solve_newsvendor(Demand, 0.1, 0.)
r2 = solve_newsvendor(Demand, 0.1, 0.5)
r3  = solve_newsvendor(Demand, 0.1, 1.)

obj = DataFrame(
a = r1[:Objective],
b = r2[:Objective],
c = r3[:Objective]
)
names!(obj, [symbol("Risk (β, λ) = (0.1, 0.0)"), symbol("Risk (β, λ) = (0.1, 0.5)") , symbol("Risk (β, λ) = (0.1, 1.0)")])

plot(
    melt(obj),
    x=:value, color=:variable,
    Geom.density,
    Theme(line_width=2pt),
    Guide.colorkey(""),
    Guide.title("max λ*E(x) + (1-λ) * CVarᵦ(x)"),
    Guide.xlabel("Profit (\$)"),
    Guide.ylabel("Frequency")
)
