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
