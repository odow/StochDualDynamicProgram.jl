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
