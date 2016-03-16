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
