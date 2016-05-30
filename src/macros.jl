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
    symin = x0.args[1]                  # name of the statein variable
    rhs   = x0.args[2]                  # values for the statein variable
    if Base.Meta.isexpr(x, :comparison) # if its a comparison
        if length(x.args) == 5          # double sided
            xin = copy(x.args[3])       # variable is in middle
        elseif length(x.args) == 3      # single comparison
            xin = copy(x.args[1])       # variable is on left
        else
            error("Unknown format for $(x)")
        end
    else
        xin = copy(x)                   # no bounds
    end
    if isa(xin, Expr)
        xin.args[1] = symin
    else
        xin = symin
    end
    code = quote end
    # create jump variables and initialise with starting solution
    push!(code.args, Expr(:(=), :stateout, Expr(:macrocall, symbol("@variable"), sp, esc(x), Expr(:kw, :start, esc(rhs)))))
    # create unbounded statein variables
    push!(code.args, Expr(:(=), :statein,  Expr(:macrocall, symbol("@variable"), sp, esc(xin))))
    push!(code.args, quote
        adddummyconstraints($sp, statein, stateout)
        stateout, statein
    end)
    return code
end

function adddummyconstraints(sp::Model, xin::JuMP.Variable, xout::JuMP.Variable)
    push!(stagedata(sp).state_vars, xout)
    push!(stagedata(sp).dual_constraints,
        @constraint(sp, xin == getvalue(xout))
    )
end
function adddummyconstraints(sp::Model, xin::Array{JuMP.Variable}, xout::Array{JuMP.Variable})
    @assert length(xin[:]) == length(xout[:])
    for i=1:length(xin[:])
        push!(stagedata(sp).state_vars, xout[i])
        push!(stagedata(sp).dual_constraints, @constraint(sp, xin[i] == getvalue(xout)[i]))
    end
end
function adddummyconstraints(sp::Model, xin::JuMP.JuMPArray, xout::JuMP.JuMPArray)
    @assert length(keys(xin)) == length(keys(xout))
    for key in keys(xin)
        push!(stagedata(sp).state_vars, xout[key...])
        push!(stagedata(sp).dual_constraints, @constraint(sp, xin[key...] == getvalue(xout[key...])))
    end
end

"""
    @scenarioconstraint(sp, [name,] rhs, constraint)

Add a scenario constraint (changes in RHS vector) to the subproblem `sp`.

Arguments:

    sp             the subproblem
    name           optional name for the scenarioconstraint
    rhs            keyword argument `key=value` where `value` is a one-dimensional array containing the scenario realisations
    constraint     any valid JuMP `@constraint` syntax that includes the keyword defined by `rhs`

Usage:

    @scenarioconstraint(sp, i=1:2, x + y <= i )
    @scenarioconstraint(sp, i=1:2, x + y <= 3 * rand(2)[i] )
    @scenarioconstraint(sp, mysc1, i=1:2, x + y <= 3 * rand(2)[i] )
"""
macro scenarioconstraint(m, args...)
    if length(args) == 3
        name = args[1]
        kw = args[2]
        c = args[3]
    elseif length(args) == 2
        name = :nothing
        kw = args[1]
        c = args[2]
    else
        error("Wrong number of arguments in @scenariocostraint")
    end

    m = esc(m)
    v = esc(kw.args[2])

    @assert length(c.args) == 3
    if c.args[2] == :(<=) || c.args[2] == :(==)
        ex = :($(c.args[1]) - $(c.args[3]))
    elseif c.args[2] == :(>=)
        ex = :($(c.args[3]) - $(c.args[1]))
    else
        error("Error in @scenarioconstraint with $c")
    end
    quote
        rhs = Float64[]
        for val in $v
            $(esc(kw.args[1])) = val
            push!(rhs, -@expression($m, $(esc(gensym())), $(esc(ex))).constant)
         end

        $(esc(kw.args[1])) = $v[1]
        con = @constraint($m, $(esc(c)))
        push!(stagedata($m).scenario_constraints, (con, rhs))
        if $(Expr(:quote, name)) != :nothing
            stagedata($m).scenario_constraint_names[$(Expr(:quote, name))] = length(stagedata($m).scenario_constraints)
        end
    end
end

"""
    @scenarioconstraints(sp, rhs, begin
        [name, ] constraint
    end)

The plural form of `@scenarioconstraint` similar to the JuMP macro `@constraints`.

Usage:

    @scenarioconstraints(sp, i=1:2, begin
               x + y <= i
               x + y <= 3 * rand(2)[i]
        mysc1, x + y <= 3 * rand(2)[i]
    end)
"""
macro scenarioconstraints(m, kw, c)
    @assert c.head == :block || error("Invalid syntax for @scenarioconstraints")
    code = quote end
    for it in c.args
        if isexpr(it, :line)
            # do nothing
        else
            if it.head == :comparison
                push!(code.args, quote
                    @scenarioconstraint($(esc(m)), $(esc(kw)), $(esc(it)))
                end)
            elseif it.head == :tuple
                if length(it.args) != 2
                    error("Unknown arguments in @scenarioconstraint")
                end
                push!(code.args, quote
                    @scenarioconstraint($(esc(m)), $(esc(it.args[1])), $(esc(kw)), $(esc(it.args[2])))
                end)
            else
                error("Unknown arguments in @scenarioconstraints")
            end
        end
    end
    push!(code.args, :(nothing))
    return code
end


"""
    @stageprofit(sp, ex)

Define the stage profit for subproblem `sp`.

Arguments:

    sp    the subproblem
    ex    a JuMP expression for the costs accrued in the subproblem

Usage:

    @stageprofit(sp, x + y)
"""
macro stageprofit(m, ex)
    m = esc(m)
    quote
        stagedata($m).stage_profit = @expression($m, $(esc(gensym())), $(esc(ex)))
    end
end
