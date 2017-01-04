#  Copyright 2017, Oscar Dowson

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
is_comparison(x) = Base.Meta.isexpr(x, :comparison) || (Base.Meta.isexpr(x, :call) && x.args[1] in comparison_symbols)

macro state(sp, x, x0)
    sp = esc(sp)                        # escape the model
    @assert x0.head == :kw              # must be a keyword
    symin, rhs = x0.args                # name of the statein variable
    if is_comparison(x)
        if length(x.args) == 5          # double sided
            xin = copy(x.args[3])       # variable is in middle
        elseif length(x.args) == 3      # single comparison
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

registerscenario!{N}(scenarios::Scenarios{N}, scenario::Scenario{N}) = push!(scenarios, scenario)
registerscenario!(sp::Model, con, values::Vector{Float64}) = registerscenario!(ext(sp).scenarios, Scenario(con, values))
registerscenario!{M,N}(scenarios::Scenarios{N}, scenario::Scenario{M}) = error("You must define scenarios with the number of options given in the SDDPModel() constructor (i.e. $(M)). You defined a scenario constraint with $(N).")

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
                    :macrocall, Symbol("@expression"),
                    sp,                         # model is first argument
                    esc(gensym()),              # generate a random Symbol
                    esc($(con.args[2]) - $(con.args[3])) # the constrexpr
                )).constant                     # want the constant term
            )
         end
        $(esc(noise_kw.args[1])) = $noise_values[1] # initialise with first scenario
        con = $(Expr(                           # add the constraint
                :macrocall, Symbol("@constraint"),
                sp,                             # the subproblem
                esc(con)                          # the constraint expression
                ))
        registerscenario!($sp, con, rhs)
        con
    end
end

# """
#     @scenario(sp, c[j=1:2], i=1:3, x <= i)
#
#     c = @constraint(sp, c[j=1:2], x <= i)
# """
# macro scenario(sp, blockname, noise, con)
# end

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
        weights        = rand(10), # probability
        (p, w) -> (p + w),         # dynamics
        (p) -> (p * x)             # Objective can use p0 as a parameter, x as a variable
    )
"""
pricestate!(sp;
    discretisation::Vector{Float64}=Float64[],
    initial::Float64=0.0,
    noise::Vector{Float64}=Float64[],
    weights::Vector{Float64}=Float64[],
    dynamics::Function=()->(),
    obj::Function=()->()
    ) = PriceScenarios(
        discretisation,
        noise,
        weights,
        dynamics,
        objective,
        0.
    )
