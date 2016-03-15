type ScenarioSet
    v::Vector{Float64}
    p::Vector{Float64}

    function ScenarioSet(Ω::Vector{Float64}, P::Vector{Float64})
        @assert sum(P) == 1 # Probabilities
        return new(Ω, P)
    end
end

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

macro addScenarioConstraint(m, kw, c)
    quote
        $(kw.args[1]) = 0
        con = @addConstraint($m, $c)
        @assert getRHS(con) == 0
        push!(m.ext[:Scenarios], (con, $(kw.args[2])))
    end
end


m=Model()
@defVar(m, x)
Ω = rand(3)
m.ext[:Scenarios] = Tuple{Any, Vector{Any}}[]
m.ext[:LastScenario] = 0
@addScenarioConstraint(m,rhs=Ω, 2x<=rhs)


display(m.linconstr[1])


function load_scenario!(sp::Model, scenario::Int)
    sp.ext[:LastScenario] = sp.ext[:CurrentScenario]
    for (c, Ω) in m.ext[:Scenarios]
        if m.ext[:LastScenario] == 0
            old_scenario = 0.
        else
            old_scenario = Ω[m.ext[:LastScenario]]
        end
        chgConstrRHS(c, getRHS(c) - old_scenario + Ω[scenario])
    end
    m.ext[:LastScenario] = scenario
    return
end

function getRHS(c::ConstraintRef{LinearConstraint})
    constr = c.m.linconstr[c.idx]
    sen = JuMP.sense(constr)
    if sen==:range
        error("Range constraints not supported for Scenarios")
    elseif sen == :(==)
        @assert constr.lb == constr.ub
        return constr.ub
    elseif sen == :>=
        return constr.lb
    elseif sen == :<=
        return constr.ub
    end
    error("Sense $(sen) not supported")
end
