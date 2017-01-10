#  Copyright 2016, Oscar Dowson
using StochDualDynamicProgram, JuMP

const sampled_errors = [-0.1290, -0.1010, -0.0814, -0.0661, -0.0530, -0.0412, -0.0303, -0.0199, -0.00987, 0.0, 0.00987, 0.0199, 0.0303, 0.0412, 0.0530, 0.0661, 0.0814, 0.1010, 0.1290]
const σ² = linspace(1, 0, 28) # Decreasing variance in changes in price over time
const transaction_cost = 0.01

const markov_states = [0.9, 1.1]
const markov_transition = [0.75 0.25; 0.3 0.7]

box(x, a, b) = min(b, max(a, x))
price_dynamics(p, w, t) = box(exp(log(p) + σ²[t]*w), 3, 9)

m = SDDPModel(
    stages     = 28,
    sense      = :Max,
    scenarios  = 19,
    transition = markov_transition
                ) do sp, t, i

    @state(sp, 0 <= contracts  <= 1.5, contracts0 = 0)
    @state(sp, 0 <= production <= 1.5, production0 = 0)

    @variables(sp, begin
        0 <= buy <= 1.2
        0 <= sell <= 1.2
        output >= 0
    end)

    @constraints(sp, begin
        contracts == contracts0 + buy - sell
        production == production0 + output
    end)

    @scenario(sp, alpha=sampled_errors, output <= alpha * markov_states[i])

    if t < 28
        price_ribs = linspace(3, 9, 5)
        price_obje = price -> (buy * price - transaction_cost * (buy + sell))
    else
        price_ribs = linspace(-5, 15, 10)
        price_obje = price -> (production - contracts) * price
    end

    objectivescenario!(sp,
        rib_locations = price_ribs,
        objective     = price_obje,
        dynamics      = (p, w) -> price_dynamics(p, w, t),
        noises        = DiscreteDistribution(
                            linspace(0, 0.05, 5),
                            fill(0.2, 5)
                        )
    )

end

# @time solvestatus = solve(m,
#     value_to_go_bound  = 1500,
#     maximum_iterations = 50,
#     policy_estimation  = MonteCarloEstimator(
#                             frequency = 1,
#                             min       = 5,
#                             max       = 100,
#                             step      = 10
#     )
#     bound_convergence  = BoundConvergence(
#                             after = 5,
#                             tol   = 1e-10
#                         ),
#     forward_pass       = ForwardPass(
#                             scenarios       = 1:10,
#                             uniformsampling = true
#                         ),
#     cut_selection      = DeMatos(5),
#     print_level        = 1
# )
#
# results = simulate(m,   # Simulate the policy
#     1000,               # number of monte carlo realisations
#     [:reservoir,        # variables to return
#     :dispatch,
#     :outflow,
#     :spill]
#     )
