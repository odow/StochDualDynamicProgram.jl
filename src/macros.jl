# This function contains JuMP extension macros
"""
Define a new state variable in the stage problem.

Usage:

    @defStateVar(m, 0<=x<=1, x0==0.5)
    @defStateVar(m, y<=1, y0==0.5)
    @defStateVar(m, z, z0==0.5)

Currently only able to handle single variables. Will break easily.
"""
macro defStateVar(m, x, x0)
    m = esc(m)
    @assert x0.head == :comparison
    if typeof(x) == Symbol
        x_sym = x
    elseif x.head == :comparison
        if length(x.args)  == 5 # Two sided
            x_sym = x.args[3]
        elseif length(x.args) == 3 # One sided
            x_sym = x.args[1]
        else
            error("Too many arguments for $(x.args)")
        end
    end
    x = esc(x)
    k = x0.args[1]
    quote
        @assert is_sp($m)
        @defVar $m $x
        @defVar $m $(esc(k))
        push!($m.ext[:state_vars], $(Expr(:quote, x_sym)))
        $m.ext[:duals][$(Expr(:quote, x_sym))] = (@addConstraint $m $(esc(x0)))
    end
end

"""
Define the value to go variable.

Usage:

    @defValueToGo(m, 0<=theta<=1)
    @defStateVar(m, theta<=1)
    @defStateVar(m, theta)

"""
macro defValueToGo(m, x)
    m = esc(m)
    if typeof(x) == Symbol
        x_sym = x
    elseif x.head == :comparison
        if length(x.args)  == 5 # Two sided
            x_sym = x.args[3]
        elseif length(x.args) == 3 # One sided
            x_sym = x.args[1]
        else
            error("Too many arguments for $(x.args)")
        end
    end
    quote
        @assert is_sp($m)
        $m.ext[:theta] = @defVar $m $(esc(x))
    end
end

"""
Right now you can only add an additive RHS with no coefficient. ie
    @addScenarioConstraint(m, rhs=[1,2,3], (...) <= (...) + rhs
"""
macro addScenarioConstraint(m, kw, c)
    m = esc(m)
    v = esc(kw.args[2])
    quote
        @assert length(collect($v)) == length($m.ext[:LastObjectives])
        $(esc(kw.args[1])) = 0
        con = @addConstraint($m, $(esc(c)))
        push!($m.ext[:Scenarios], (con, collect($v)))
    end
end

macro setStageProfit(m, ex)
    m = esc(m)
    quote
        $m.ext[:StageProfit] = @defExpr $(esc(gensym())) $(esc(ex))
    end
end
