# Copyright 2017, Oscar Dowson

setvalue!(x::StateVariable, y::Float64) = JuMP.setrhs!(x.c, y)
JuMP.getdual(x::StateVariable) = JuMP.getdual(x.c)

function statevariable!(states::Vector{StateVariable}, xin::JuMP.Variable, xout::JuMP.Variable)
    push!(states,
            StateVariable(
                xout,
                @constraint(sp, xin == getvalue(xout))
            )
    )
end

function statevariable!(states::Vector{StateVariable}, xin::Array{JuMP.Variable}, xout::Array{JuMP.Variable})
    @assert length(xin) == length(xout)
    for i in 1:length(xin)
        statevariable!(states, xin[i], xout[i])
    end
end

function statevariable!{T<:Union{JuMP.JuMPArray, JuMP.JuMPDict}}(states::Vector{StateVariable}, xin::T, xout::T)
    @assert length(keys(xin)) == length(keys(xout))
    for key in keys(xin)
        statevariable!(states, xin[key...], xout[key...])
    end
end

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
macro state(sp, x, x0)
    sp = esc(sp)                        # escape the model
    @assert x0.head == :kw              # must be a keyword
    symin, rhs = x0.args                # name of the statein variable
    if is_comparison(x)
        if length(x.args) == 5          # double sided
            xin = identity(x.args[3])       # variable is in middle
        elseif length(x.args) == 3      # single comparison
            xin = identity(x.args[2])       # variable is second entry
        else
            error("Unknown format for $(x)")
        end
    else
        xin = identity(x)                   # no bounds
    end
    if isa(xin, Expr)                   # x has indices
        xin.args[1] = symin             # so just change the name
    else                                # its just a Symbol
        xin = symin                     # so change the Symbol
    end
    quote
        stateout = $(Expr(:macrocall, Symbol("@variable"), sp, esc(x), Expr(:kw, :start, esc(rhs))))
        statein  = $(Expr(:macrocall, Symbol("@variable"), sp, esc(xin)))
        statevariable!(states($sp), statein, stateout)
        stateout, statein
    end
end
