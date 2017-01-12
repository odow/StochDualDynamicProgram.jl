# Copyright 2017, Oscar Dowson

function scenario!(scenarios::Vector{Scenario}, con::LinearConstraint, values::Vector{Float64})
    @assert length(scenarios) == length(values)
    for i in 1:length(values)
        push!(scenarios[i].arr, ConstraintRHS(con, values[i]))
    end
end

"""
    @scenario(sp, i=1:3, x <= i)
"""
macro scenario(sp, noise_kw, con)
    sp = esc(sp)
    noise_values = esc(noise_kw.args[2])
    quote
        rhs = Float64[]                         # intialise RHS vector
        for noise in $noise_values    # for each scenario
            $(esc(noise_kw.args[1])) = noise  # set the scenariovalue
            push!(rhs,                          # add to the rhs vector
                -$(Expr(                        # negate to shift from LHS to RHS
                    :call, :(-),
                    esc(con.args[2]),                         # model is first argument
                    esc(con.args[3]) # the constrexpr
                )).constant                     # want the constant term
            )
         end
        $(esc(noise_kw.args[1])) = $noise_values[1] # initialise with first scenario
        full_con = $(Expr(                           # add the constraint
                :macrocall, Symbol("@constraint"),
                sp,                             # the subproblem
                esc(con)                          # the constraint expression
                ))
        scenario!(scenarios($sp), full_con, rhs)
        full_con
    end
end

"""
    @scenarios(sp, i=1:3, begin
        x <= i
        y <= A[i]
    end)
"""
macro scenarios(sp, noise, blkargs)
    @assert blkargs.head == :block || error("Invalid syntax for @scenarioconstraints")
    code = quote end
    for con in blkargs.args
        if !Base.Meta.isexpr(con, :line)
            if is_comparison(con)
                push!(code.args, Expr(:macrocall, Symbol("@scenario"), esc(sp), esc(noise), esc(con)))
            else
                error("Unknown arguments in @scenarios")
            end
        end
    end
    push!(code.args, :(nothing))
    return code
end
