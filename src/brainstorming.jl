# macro addScenarioConstraint(m, kw, c)
#     quote
#         $(kw.args[1]) = 1
#         @defExpr(ex1, $(c.args[1]) - $(c.args[3]))
#         $(kw.args[1]) = 0
#         @defExpr(ex0, $(c.args[1]) - $(c.args[3]))
#
#         con = @addConstraint($m, $c)
#         push!(m.ext[:scenario_constraints], (con, -(ex1 - ex0).constant, $(kw.args[2])))
#     end
# end
m=Model()
@defVar(m, x)
Ω = rand(3)
m.ext[:scenario_constraints] = Tuple{Any, Vector{Any}}[]
m.ext[:LastScenario] = 0
@addScenarioConstraint(m,rhs=Ω, 2x<=rhs)
display(m.linconstr[1])


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
        $m.ext[:dual_constraints][$(Expr(:quote, x_sym))] = (@addConstraint $m $(esc(x0)))
    end
end

workspace()
using JuMP

macro foobar(x)
    @show x
    @show x.args
    @show x.args[2]
    @show x.args[2].args
    return
end
@foobar X[r=1:3]

macro foo1(m, x, x0)
    simple_var = true
    m = esc(m)
    if typeof(x) == Symbol
        x_sym = x
    elseif typeof(x) == Expr
        if x.head == :comparison
            if length(x.args)  == 5 # Two sided
                x_sym = x.args[3]
            elseif length(x.args) == 3 # One sided
                x_sym = x.args[1]
            else
                error("Too many arguments for $(x.args)")
            end
        elseif x.head == :ref
            x_sym = x.args[1]
            simple_var = false
        else
            error("$(x) not recognised")
        end
    end
    x = esc(x)

    @assert x0.head == :comparison
    @assert length(x0.args) == 3
    kl = esc(x0.args[1])
    kr = esc(x0.args[3])
    if simple_var
        quote
            # @assert is_sp($m)
            @defVar $m $x
            @defVar $m $kl
            push!($m.ext[:state_vars], $(Expr(:quote, x_sym)))
            $m.ext[:dual_constraints][$(Expr(:quote, x_sym))] = (@addConstraint $m $(esc(x0)))
        end
    else
        quote
            # @assert is_sp($m)
            a = @defVar $m $x
            push!($m.ext[:state_vars], ($(Expr(:quote, x_sym)), keys(a)))
            b = @defVar $m $kl
            @show b
            @show $(esc(x0.args))
            @assert length(keys(a)) == length(intersect(keys(a), keys(b)))
            for key in keys(a)
                $m.ext[:dual_constraints][($(Expr(:quote, x_sym)), key)] = @addConstraint($m, $(esc(x0)))#b[key...] == $kr[key...])
            end
        end
    end
end

m = Model(); m.ext[:state_vars] = Any[]; m.ext[:dual_constraints] = Dict{Any, Any}();T = [:upper, :lower]; R = Dict{Symbol, Float64}(:upper=>1., :lower=>2.); @foo1(m, y[t=T], y0[t=T]==R[t])
m.ext[:state_vars]
m.ext[:dual_constraints]

m = Model(); m.ext[:state_vars] = Any[]; m.ext[:dual_constraints] = Dict{Any, Any}();@foo1(m, y, y0==1)
m.ext[:state_vars]
m.ext[:dual_constraints]


function defStateVariable!{T}(m::Model, I::Vector{T})#, lb::Any, ub::Any, initial::Any)
    s = gensym()
    push!(m.ext[:state_vars], s)
    x = @defVar(m, 0 <= s[i=collect([:upper, :lower])] <= 1.)
    s0 = gensym()
    x0 = @defVar(m, s0)
    # m.ext[:dual_constraints][s] = @addConstraint(m, x0 == 1)#initial[i])
    return x, x0
end
m=Model();m.ext[:state_vars] = Any[];m.ext[:dual_constraints] = Any[]
x, x0 = defStateVariable!(m, [:upper, :lower])

x

function StateVariable(m::Model, lb::Number, ub::Number, initial::Number, name::String)
        s = gensym()
        push!(m.ext[:state_vars], s)
        x = @defVar(m, lb[i] <= s <= ub[i])
        s0 = gensym()
        x0 = @defVar(m, s0)
        m.ext[:dual_constraints][s] = @addConstraint(m, x0 == initial[i])
    return x, x0
end

reservoir_max = Dict{Symbol, Float64}(
    :upper => 200,
    :lower => 200
    )


    reservoir_max = Dict{Symbol, Float64}(
        :upper => 200,
        :lower => 200
        )

# Initial fill
reservoir_initial = Dict{Symbol, Float64}(
    :upper => 200,
    :lower => 200
    )
x, x0 = defStateVariable!(m, [:upper, :lower], [0], reservoir_max, reservoir_initial)
















println("The end.")
