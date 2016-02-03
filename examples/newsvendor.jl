using StochDualDynamicProgram, JuMP

# Hydro problem with different transition matrix
function solve_newsvendor(;
    Demand = [
        10 15;
        12 20;
        8  20
    ],

    PurchasePrice = [5, 8],

    RetailPrice = 7,

    # Transition matrix
    Transition = Array{Float64, 2}[
        [0.6 0.4; 0.3 0.7],
        [0.3 0.7; 0.3 0.7],
        [0.5 0.5; 0.5 0.5]
      ]

      )

    # Initialise SDDP Model
    m = SDDPModel(stages=3, markov_states=2, transition=Transition, initial_markov_state=1) do sp, stage, markov_state
        @defStateVar(sp, 0 <= stock <= 100, stock0==5)

        @defVar(sp, buy>=0)
        @defVar(sp, 0 <= sell <= Demand[stage, markov_state])

        if stage < 3
            @defValueToGo(sp, theta <= 1000)
            @setObjective(sp, Max, theta + sell * RetailPrice - buy * PurchasePrice[markov_state])
        else
            @setObjective(sp, Max, sell * RetailPrice - buy * PurchasePrice[markov_state])
        end

        @addConstraint(sp, stock == stock0 + buy - sell)
    end

    solve(m,                # Solve the model using the SDDP algorithm
        forward_passes=1000,  # number of realisations in bound simulation
        backward_passes=10   # number of cutting iterations before convergence check
    )

    results = simulate(m,   # Simulate the policy
        1000,               # number of monte carlo realisations
        [:stock, :buy, :sell]
        )

    return results
end
