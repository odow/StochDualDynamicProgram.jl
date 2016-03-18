# macro addScenarioConstraint(m, kw, c)
#     quote
#         $(kw.args[1]) = 1
#         @defExpr(ex1, $(c.args[1]) - $(c.args[3]))
#         $(kw.args[1]) = 0
#         @defExpr(ex0, $(c.args[1]) - $(c.args[3]))
#
#         con = @addConstraint($m, $c)
#         push!(m.ext[:Scenarios], (con, -(ex1 - ex0).constant, $(kw.args[2])))
#     end
# end
m=Model()
@defVar(m, x)
Ω = rand(3)
m.ext[:Scenarios] = Tuple{Any, Vector{Any}}[]
m.ext[:LastScenario] = 0
@addScenarioConstraint(m,rhs=Ω, 2x<=rhs)
display(m.linconstr[1])

# Thinking about risk aversion

function risk_averse_weightings(x::Vector{Float64}, p::Vector{Float64},  ß::Float64=0.5, ismax::Bool=true)
    n = 10
    I = sortperm(x, rev=!ismax)
    y = similar(x)
    q = 0.
    for i in I
        q >=  ß && break
        y[i] = min(p[i],  ß - q)
        q += y[i]
    end
    return y ./ ß
end
risk_averse_weightings(x::Vector{Float64}, beta::Float64=0.5) = risk_averse_weightings(x, ones(length(x)) / length(x), beta)


x = rand(10)
y = risk_averse_weightings(x, 0.5)
dot(x, y)
mean(x)


function add_cut!{M,N,S,T}(m::SDDPModel{M,N,S,T}, stage::Int, markov::Int)
    sp = m.stage_problems[stage, markov]

    P = zeros(S*M)
    i=1
    for mkv=1:M
        for s=1:S
            P[i] = get_transition(m, stage, markov, mkv)/S
            i+=1
        end
    end
    w = risk_averse_weightings([sp.ext[:LastObjectives] for sp in m.stage_problems[stage+1, :]], P,  m.beta_quantile)

    Prob = zeros(M,S)
    i=1
    for mkv=1:M
        for s=1:S
            Prob[mkv, s] = P[i] * w[i]
            i+=1
        end
    end

    @defExpr(rhs, sum{
        sum{
            Prob[mkv, s] * (
                m.stage_problems[stage+1, mkv].ext[:LastObjectives][s] +
                sum{
                    m.stage_problems[stage+1, mkv].ext[:DualValues][state][s] * (
                        getVar(sp, state) - getRHS(m.stage_problems[stage+1, mkv].ext[:duals][state])
                    )
                , state in sp.ext[:state_vars]}
            )
        ,s=1:S}
    , mkv=1:N}
    )

    if m.sense==:Max
        @addConstraint(sp, sp.ext[:theta] <= rhs)
    else
        @addConstraint(sp, sp.ext[:theta] >= rhs)
    end
end






function set_valid_bound!{M,N,S,T}(m::SDDPModel{M,N,S,T})
    obj = 0.0
    if m.init_markov_state == 0

        P = zeros(S*N)
        i=1
        for mkv=1:N
            for s=1:S
                P[i] = get_transition(m, 1, markov, mkv)/S
                i+=1
            end
        end
        w = risk_averse_weightings(vcat([sp.ext[:LastObjectives] for sp in m.stage_problems[1, :]]...), P,  m.beta_quantile, m.sense==:Max)
        i=1
        for mkv=1:N
            for s=1:S
                obj += w[i] * m.stage_problems[1, mkv].ext[:LastObjectives][s]
                i+=1
            end
        end

    else

        P = ones(S) / S
        w = risk_averse_weightings(m.stage_problems[1, m.init_markov_state]sp.ext[:LastObjectives], P,  m.beta_quantile, m.sense==:Max)
        for s=1:S
            obj += w[s] * m.stage_problems[1, m.init_markov_state].ext[:LastObjectives][s]
        end

    end


    if m.sense==:Max && obj < getBound(m)
        setBound!(m, obj)
    elseif m.sense==:Min && obj > getBound(m)
        setBound!(m, obj)
    end
end



function t_test(x::Vector; conf_level=0.95)
    tstar = quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    SE = std(x)/sqrt(length(x))
    lo, hi = mean(x) + [-1, 1] * tstar * SE
    return (lo, hi)#, mean(x))
end

using Distributions

function cvar_estimate(x, beta, lambda, n=1000)
    y = zeros(100)#div(length(x), 100))
    z = zeros(n)
    for i=1:n
        rand!(y, x)
        z[i] = cvar(y, beta, lambda)
    end
    t_test(z)
end

function cvar{T}(x::Vector{T}, beta::Float64=1., lambda::Float64=1.)
    @assert beta >= 0 && beta <= 1.
    @assert lambda >= 0 && lambda <= 1.
    lambda * mean(x) + (1 - lambda) * mean(x[x.<quantile(x, beta)])
end

using Gadfly

x = randn(2000)
b = 0.05
l = 0.


function plotstuff(n)
    plot(
    layer(xintercept=[cvar(x, b, l)], Geom.vline, Theme(line_width=2pt, default_color=colorant"slategray")),
    layer(x=cvar_estimate(x, b, l,n), Geom.histogram, Theme(default_color=colorant"red")),
    layer(x=x, Geom.histogram)
    )
end
plotstuff(1000)















using StochDualDynamicProgram, JuMP

# Hydro problem with different transition matrix
function solve_newsvendor(Demand, beta_quant=0.5)
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
            theta_bound = 1000
        ) do sp, stage, markov_state

        # ====================
        #   State variable
        @defStateVar(sp, 0 <= stock <= 100, stock0==5)

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
        forward_passes=1000,  # number of realisations in bound simulation
        backward_passes=10,  # number of cutting iterations before convergence check
        beta_quantile=beta_quant,
        risk_lambda = 0.5,
        max_iters=100
    )

    results = simulate(m,   # Simulate the policy
        3000,               # number of monte carlo realisations
        [:stock, :buy, :sell]
        )

    return results
end


Demand = 20 * rand(3,30)

r_risk = solve_newsvendor(Demand, 0.1);
r_exp  = solve_newsvendor(Demand, 1.);

plot(
layer(x=r_risk[:Objective], Geom.density, Theme(default_color=colorant"red")),
layer(x=r_exp[:Objective], Geom.density)
)
