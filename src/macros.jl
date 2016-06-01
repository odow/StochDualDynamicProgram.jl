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
    if isa(xin, Expr)                   # x has indices
        xin.args[1] = symin             # so just change the name
    else                                # its just a symbol
        xin = symin                     # so change the symbol
    end
    quote
        stateout = $(Expr(:macrocall, symbol("@variable"), sp, esc(x), Expr(:kw, :start, esc(rhs))))
        statein  = $(Expr(:macrocall, symbol("@variable"), sp, esc(xin)))
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
    @assert c.head == :comparison               # check c is a comparison constraint
    @assert length(c.args) == 3                 # check that it has (LHS, (comparison), RHS)
    @assert c.args[2]  in [:(<=), :(==), :(>=)] # check valid constraint type
    constrexpr = :($(c.args[1]) - $(c.args[3])) # LHS - RHS
    quote
        rhs = Float64[]                         # intialise RHS vector
        for scenariovalue in $scenariovalues    # for each scenario
            $(esc(kw.args[1])) = scenariovalue  # set the scenariovalue
            push!(rhs,                          # add to the rhs vector
                -$(Expr(                        # negate to shift from LHS to RHS
                    :macrocall, symbol("@expression"),
                    sp,                         # model is first argument
                    esc(gensym()),              # generate a random symbol
                    esc(constrexpr)             # the constrexpr
                )).constant                     # want the constant term
            )
         end

        $(esc(kw.args[1])) = $scenariovalues[1] # initialise with first scenario
        con = $(Expr(                           # add the constraint
                :macrocall, symbol("@constraint"),
                sp,                             # the subproblem
                esc(c)                          # the constraint expression
                ))
        registerscenarioconstraint!($sp, con, rhs, $(Expr(:quote, name)))
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
            if it.head == :comparison
                push!(code.args,
                    Expr(:macrocall, symbol("@scenarioconstraint"), esc(m), esc(kw), esc(it))
                )
            elseif it.head == :tuple
                if length(it.args) != 2
                    error("Unknown arguments in @scenarioconstraint")
                end
                push!(code.args,
                    Expr(:macrocall, symbol("@scenarioconstraint"), esc(m), esc(it.args[1]), esc(kw), esc(it.args[2]))
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
