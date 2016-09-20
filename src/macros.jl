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
    push!(stagedata(sp).state_vars, xout)
    push!(stagedata(sp).dual_constraints, @constraint(sp, xin == getvalue(xout)))
end
function registerstatevariable!(sp::Model, xin::Array{JuMP.Variable}, xout::Array{JuMP.Variable})
    @assert length(xin[:]) == length(xout[:])
    for i=1:length(xin[:])
        registerstatevariable!(sp, xin[i], xout[i])
    end
end
function registerstatevariable!(sp::Model, xin::JuMP.JuMPArray, xout::JuMP.JuMPArray)
    @assert length(keys(xin)) == length(keys(xout))
    for key in keys(xin)
        registerstatevariable!(sp, xin[key...], xout[key...])
    end
end
function registerstatevariable!(sp::Model, xin::JuMP.JuMPDict, xout::JuMP.JuMPDict)
    @assert length(keys(xin)) == length(keys(xout))
    for key in keys(xin)
        registerstatevariable!(sp, xin[key...], xout[key...])
    end
end

macro states(m, b)
    @assert b.head == :block || error("Invalid syntax for @states")
    code = quote end
    for it in b.args
        if Base.Meta.isexpr(it, :line)
            # do nothing
        else
            if it.head == :tuple
                if length(it.args) != 2
                    error("Unknown arguments in @states")
                end
                if it.args[2].head == :(=)
                    it.args[2].head = :kw
                end
                push!(code.args,
                    Expr(:macrocall, Symbol("@state"), esc(m), esc(it.args[1]), esc(it.args[2]))
                )
            else
                error("Unknown arguments in @states")
            end
        end
    end
    push!(code.args, :(nothing))
    return code
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
macro scenarioconstraint(sp, args...)
    sp = esc(sp)                                # escape the model
    if length(args) == 3                        # the named version
        name, kw, c = args                      # unpack the arguments
    elseif length(args) == 2                    # no name given
        name = :nothing                         # blank name
        kw, c = args                            # unpack the arguments
    else                                        # need 2 or 3 arguemnts
        error("Wrong number of arguments in @scenariocostraint")
    end
    @assert kw.head == :kw                      # check its a keyword
    scenariovalues = esc(kw.args[2])            # get the vector of values
    @assert c.head == :call               # check c is a comparison constraint
    @assert length(c.args) == 3                 # check that it has (LHS, (comparison), RHS)
    @assert c.args[1]  in comparison_symbols # check valid constraint type
    constrexpr = :($(c.args[2]) - $(c.args[3])) # LHS - RHS
    quote
        rhs = Float64[]                         # intialise RHS vector
        for scenariovalue in $scenariovalues    # for each scenario
            $(esc(kw.args[1])) = scenariovalue  # set the scenariovalue
            push!(rhs,                          # add to the rhs vector
                -$(Expr(                        # negate to shift from LHS to RHS
                    :macrocall, Symbol("@expression"),
                    sp,                         # model is first argument
                    esc(gensym()),              # generate a random Symbol
                    esc(constrexpr)             # the constrexpr
                )).constant                     # want the constant term
            )
         end

        $(esc(kw.args[1])) = $scenariovalues[1] # initialise with first scenario
        con = $(Expr(                           # add the constraint
                :macrocall, Symbol("@constraint"),
                sp,                             # the subproblem
                esc(c)                          # the constraint expression
                ))
        registerscenarioconstraint!($sp, con, rhs, $(Expr(:quote, name)))
        con
    end
end

function registerscenarioconstraint!(sp::Model, con, rhs, name)
    push!(stagedata(sp).scenario_constraints, (con, rhs))
    if name != :nothing
        stagedata(sp).scenario_constraint_names[name] = length(stagedata(sp).scenario_constraints)
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
        if Base.Meta.isexpr(it, :line)
            # do nothing
        else
            if it.head == :call && it.args[1] in comparison_symbols
                push!(code.args,
                    Expr(:macrocall, Symbol("@scenarioconstraint"), esc(m), esc(kw), esc(it))
                )
            elseif it.head == :tuple
                if length(it.args) != 2
                    error("Unknown arguments in @scenarioconstraint")
                end
                push!(code.args,
                    Expr(:macrocall, Symbol("@scenarioconstraint"), esc(m), esc(it.args[1]), esc(kw), esc(it.args[2]))
                )
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
function addpenalties!(sp, weight)
    @variable(sp, _penalty_upper_ >= 0)
    @variable(sp, _penalty_lower_ >= 0)
    stagedata(sp).relaxationpenalties += weight * (_penalty_lower_ + _penalty_upper_)
    return _penalty_upper_, _penalty_lower_
end
macro relaxedconstraint(sp, weight, args...)
    code = quote end
    sp = esc(sp)
    if length(args) == 1
        constr = args[1]
    elseif length(args) == 2
        constr = args[2]
    else
        error("Too many arguments in @relaxedconstraint")
    end
    @assert constr.head == :call && constr.args[1] in comparison_symbols
    _penalty_upper_, _penalty_lower_ = gensym(), gensym()
    push!(code.args, quote
        $(esc(_penalty_upper_)), $(esc(_penalty_lower_)) = addpenalties!($sp, $weight)
    end)
    if constr.args[1] == :(<=)
        newconstr = Expr(:call, :(<=), Expr(:call, :(-), constr.args[2], _penalty_lower_), constr.args[3])
    elseif constr.args[1] == :(>=)
         newconstr = Expr(:call, :(>=), Expr(:call, :(+), constr.args[2], _penalty_upper_), constr.args[3])
    elseif constr.args[1] == :(==)
         newconstr = Expr(:call, :(==), Expr(:call, :(-), Expr(:call, :(+), constr.args[2], _penalty_upper_), _penalty_lower_), constr.args[3])
    else
        error("Comparison operator $(constr) not recognised")
    end
    if length(args) == 1
        push!(code.args, Expr(:macrocall, Symbol("@constraint"), sp, esc(newconstr)))
    else
        push!(code.args,Expr(:macrocall, Symbol("@constraint"), sp, esc(args[1]), esc(newconstr)))
    end
    return code
end

macro relaxedconstraints(m, penalty, c)
    @assert c.head == :block || error("Invalid syntax for @relaxedconstraints")
    code = quote end
    for it in c.args
        if Base.Meta.isexpr(it, :line)
            # do nothing
        else
            if it.head == :call && it.args[1] in comparison_symbols
            # if it.head == :comparison
                push!(code.args,
                    Expr(:macrocall, Symbol("@relaxedconstraint"), esc(m), esc(penalty), esc(it))
                )
            elseif it.head == :tuple
                if length(it.args) != 2
                    error("Unknown arguments in @relaxedconstraints")
                end
                push!(code.args,
                    Expr(:macrocall, Symbol("@relaxedconstraint"), esc(m), esc(penalty), esc(it.args[1]), esc(it.args[2]))
                )
            else
                error("Unknown arguments in @relaxedconstraints")
            end
        end
    end
    push!(code.args, :(nothing))
    return code
end
