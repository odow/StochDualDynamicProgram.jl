#  Copyright 2016, Oscar Dowson

"""
    @state(sp, stateleaving, stateentering)

Define a new state variable in the subproblem `sp`.

Arguments:

    sp               the subproblem
    stateleaving     any valid JuMP `@variable` syntax to define the value of the state variable at the end of the stage
    stateentering    any valid JuMP `@variable` syntax to define the value of the state variable at the beginning of the stage

Usage:

    @state(sp, 0 <= x[i=1:3] <= 1, x0=rand(3)[i] )
    @state(sp,      y        <= 1, y0=0.5        )
    @state(sp,      z            , z0=0.5        )

"""
Base.copy(s::Symbol) = s
const comparison_symbols = [:(<=), :(>=), :(==)]
macro state(sp, x, x0)
    sp = esc(sp)                        # escape the model
    @assert x0.head == :kw              # must be a keyword
    symin, rhs = x0.args                # name of the statein variable
    if Base.Meta.isexpr(x, :comparison)
        if length(x.args) == 5          # double sided
            xin = copy(x.args[3])       # variable is in middle
        end
    elseif Base.Meta.isexpr(x, :call) && x.args[1] in comparison_symbols # if its a comparison
        if length(x.args) == 3      # single comparison
            xin = copy(x.args[2])       # variable is second entry
        else
            error("Unknown format for $(x)")
        end
    else
        xin = copy(x)                   # no bounds
    end
    if isa(xin, Expr)                   # x has indices
        xin.args[1] = symin             # so just change the name
    else                                # its just a Symbol
        xin = symin                     # so change the Symbol
    end
    quote
        stateout = $(Expr(:macrocall, Symbol("@variable"), sp, esc(x), Expr(:kw, :start, esc(rhs))))
        statein  = $(Expr(:macrocall, Symbol("@variable"), sp, esc(xin)))
        registerstatevariable!($sp, statein, stateout)
        stateout, statein
    end
end

function registerstatevariable!(sp::Model, xin::JuMP.Variable, xout::JuMP.Variable)
    push!(ext(sp).states,
            StateVariable(
                xout,
                @constraint(sp, xin == getvalue(xout))
            )
    )
end

function registerstatevariable!(sp::Model, xin::Array{JuMP.Variable}, xout::Array{JuMP.Variable})
    @assert length(xin[:]) == length(xout[:])
    for i=1:length(xin[:])
        registerstatevariable!(sp, xin[i], xout[i])
    end
end
function registerstatevariable!{T<:Union{JuMP.JuMPArray, JuMP.JuMPDict}}(sp::Model, xin::T, xout::T)
    @assert length(keys(xin)) == length(keys(xout))
    map(key->registerstatevariable!(sp, xin[key...], xout[key...]), keys(xin))
end

"""
    @scenario(sp, i=1:3, x <= i)
"""
macro scenario(sp, noise, con)
end

"""
    @scenarios(sp, i=1:3, begin
        x <= i
        y <= A[i]
    end)
"""
macro scenarios(sp, noise, args)
end

"""
    @stageobjective(sp, ex)

Define the stage profit for subproblem `sp`.

Arguments:

    sp    the subproblem
    ex    a JuMP expression for the costs accrued in the subproblem

Usage:

    @stageobjective(sp, x + y)
"""
macro stageobjective(m, ex)
    m = esc(m)
    quote
        ext($m).stageobjective = @expression($m, $(esc(gensym())), $(esc(ex)))
    end
end

"""
    pricestate!(sp,              # subproblem
        discretisation = 0:10,     # outgoing price
        initial        = 0.5,      # incoming price
        noise          = rand(10), # noises
        (p, w) -> (p + w),         # dynamics
        (p) -> (p * x)             # Objective can use p0 as a parameter, x as a variable
    )
"""
function pricestate!(sp, discretisaion, initial, noise, dynamics::Function, obj::Function)
end
