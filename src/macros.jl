# This function contains JuMP extension macros
# """
# Define a new state variable in the stage problem.
#
# Usage:
#
#     @defStateVar(m, 0<=x<=1, x0==0.5)
#     @defStateVar(m, y<=1, y0==0.5)
#     @defStateVar(m, z, z0==0.5)
#
# Currently only able to handle single variables. Will break easily.
# """

# Much of the functionality in @defStateVar is shamelessly copied from
# both the JuMP @defVar and JuMPer @defUnc(https://github.com/IainNZ/JuMPeR.jl/blob/master/src/robustmacro.jl)

import JuMP: assert_validmodel, validmodel, esc_nonconstant
import JuMP: getloopedcode, buildrefsets, getname, registervar
import JuMP: storecontainerdata, isdependent, JuMPContainerData, pushmeta!, JuMPContainer
import JuMP: EMPTYSTRING
using Base.Meta

function State(m::Model, lower::Number, upper::Number, name)
    v = Variable(m,lower,upper,:Cont,utf8(string(name)),NaN)
    push!(stagedata(m).state_vars, v)
    return v
end

function State0(m::Model, init, name, name0)
    v0 = Variable(m,-Inf,Inf,:Cont,utf8(string(name0)),init)
    push!(stagedata(m).dual_constraints, @addConstraint(m, v0 == init))
    return v0
end


macro defStateVar(args...)
    length(args) <= 1 &&
        error("in @defStateVar: expected model as first argument, then variable information.")
    m = esc(args[1])
    x = args[2]
    x0 = args[3]
    extra = vcat(args[4:end]...)

    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    if isexpr(x,:comparison)
        # We have some bounds
        if x.args[2] == :>= || x.args[2] == :≥
            if length(x.args) == 5
                # ub >= x >= lb
                x.args[4] == :>= || x.args[4] == :≥ || error("Invalid variable bounds")
                var = x.args[3]
                lb = esc_nonconstant(x.args[5])
                ub = esc_nonconstant(x.args[1])
            else
                # x >= lb
                var = x.args[1]
                @assert length(x.args) == 3
                lb = esc_nonconstant(x.args[3])
                ub = Inf
            end
        elseif x.args[2] == :<= || x.args[2] == :≤
            if length(x.args) == 5
                # lb <= x <= u
                var = x.args[3]
                (x.args[4] != :<= && x.args[4] != :≤) &&
                    error("in @defStateVar ($var): expected <= operator after variable name.")
                lb = esc_nonconstant(x.args[1])
                ub = esc_nonconstant(x.args[5])
            else
                # x <= ub
                var = x.args[1]
                # NB: May also be lb <= x, which we do not support
                #     We handle this later in the macro
                @assert length(x.args) == 3
                ub = esc_nonconstant(x.args[3])
                lb = -Inf
            end
        elseif x.args[2] == :(==)
            # fixed variable
            var = x.args[1]
            @assert length(x.args) == 3
            lb = esc(x.args[3])
            ub = esc(x.args[3])
        else
            # Its a comparsion, but not using <= ... <=
            error("in @defStateVar ($(string(x))): use the form lb <= ... <= ub.")
        end
    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        lb = -Inf
        ub = Inf
    end

    # separate out keyword arguments
    extra = filter(ex->!isexpr(ex,:kw), extra)

    # process keyword arguments
    x0_name = x0.args[1]
    x0_value = esc(x0.args[2])

    quotx0name = quot(getname(x0_name))
    escx0name = esc(getname(x0_name))

    quotvarname = quot(getname(var))
    escvarname  = esc(getname(var))

    if isa(var,Symbol)
        # Easy case - a single variable
        return assert_validmodel(m, quote
            $(esc(var)) = State($m, $lb, $ub, $quotvarname)
            registervar($m, $quotvarname, $escvarname)
            # push!($stagedata(m).state_vars, $quotvarname)

            $(esc(x0_name)) = State0($m, $x0_value, $quotvarname, $quotx0name)
            registervar($m, $quotx0name, $escx0name)

            $(esc(var)), $(esc(x0_name))
        end)
    end
    isa(var,Expr) || error("in @defStateVar: expected $var to be a variable name")

    var2 = deepcopy(var)
    var2.args[1] = x0_name

    quotx0name = quot(getname(var2))
    escx0name = esc(getname(var2))

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall1, idxvars1, idxsets1, idxpairs1, condition1 = buildrefsets(var)
    clear_dependencies1(i) = (isdependent(idxvars1,idxsets1[i],i) ? nothing : idxsets1[i])

    refcall2, idxvars2, idxsets2, idxpairs2, condition2 = buildrefsets(var2)
    clear_dependencies2(i) = (isdependent(idxvars2,idxsets2[i],i) ? nothing : idxsets2[i])
    # code = :( $(refcall) = State($m, $lb, $ub, $x0_value, EMPTYSTRING, EMPTYSTRING) )
    code1 = :( $(refcall1) = State($m, $lb, $ub, $quotvarname) )
    code2 = :( $(refcall2) = State0($m, $x0_value, $quotvarname, $quotx0name) )

    looped1 = getloopedcode(var, code1, condition1, idxvars1, idxsets1, idxpairs1, :Variable)
    looped2 = getloopedcode(var2, code2, condition2, idxvars2, idxsets2, idxpairs2, :Variable)

    return assert_validmodel(m, quote
        $looped1
        push!($(m).dictList, $escvarname)
        # push!($stagedata(m).state_vars, $quotvarname)
        registervar($m, $quotvarname, $escvarname)
        storecontainerdata($m, $escvarname, $quotvarname,
                           $(Expr(:tuple,map(clear_dependencies1,1:length(idxsets1))...)),
                           $idxpairs1, $(quot(condition1)))

        $looped2
        push!($(m).dictList, $escx0name)
        registervar($m, $quotx0name, $escx0name)
        storecontainerdata($m, $escx0name, $quotx0name,
                          $(Expr(:tuple,map(clear_dependencies2,1:length(idxsets2))...)),
                          $idxpairs2, $(quot(condition2)))

        isa($escvarname, JuMPContainer) && pushmeta!($escvarname, :model, $m)
        isa($escx0name, JuMPContainer) && pushmeta!($escx0name, :model, $m)

        $escvarname, $escx0name
    end)

end

# macro defStateVar(m, x, x0)
#     m = esc(m)
#     @assert x0.head == :comparison
#     if typeof(x) == Symbol
#         x_sym = x
#     elseif x.head == :comparison
#         if length(x.args)  == 5 # Two sided
#             x_sym = x.args[3]
#         elseif length(x.args) == 3 # One sided
#             x_sym = x.args[1]
#         else
#             error("Too many arguments for $(x.args)")
#         end
#     end
#     x = esc(x)
#     k = x0.args[1]
#     quote
#         @assert is_sp($m)
#         @defVar $m $x
#         @defVar $m $(esc(k))
#         push!($stagedata(m).state_vars, $(Expr(:quote, x_sym)))
#         $stagedata(m).dual_constraints[$(Expr(:quote, x_sym))] = (@addConstraint $m $(esc(x0)))
#     end
# end

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
        stagedata($m).theta = @defVar $m $(esc(x))
    end
end

"""
Right now you can only add an additive RHS with no coefficient. ie
    @addScenarioConstraint(m, rhs=[1,2,3], (...) <= (...) + rhs
"""
# macro addScenarioConstraint(m, kw, c)
#     m = esc(m)
#     v = esc(kw.args[2])
#     quote
#         @assert length(collect($v)) == length(stagedata($m).objective_value)
#         $(esc(kw.args[1])) = 0
#         con = @addConstraint($m, $(esc(c)))
#         push!(stagedata($m).scenario_constraints, (con, collect($v)))
#     end
# end

macro addScenarioConstraint(m, kw, c)
    m = esc(m)
    v = esc(kw.args[2])

    @assert length(c.args) == 3
    if c.args[2] == :(<=) || c.args[2] == :(==)
        ex = :($(c.args[1]) - $(c.args[3]))
    elseif c.args[2] == :(>=)
        ex = :($(c.args[3]) - $(c.args[1]))
    else
        error("Error in @addScenarioConstraint with $c")
    end
    quote
        rhs = Float64[]
        for val in $v
            $(esc(kw.args[1])) = val
            push!(rhs, -@defExpr($(esc(ex))).constant)
         end

        $(esc(kw.args[1])) = $v[1]
        con = @addConstraint($m, $(esc(c)))
        push!(stagedata($m).scenario_constraints, (con, rhs))
    end
end
# @addScenarioConstraint2(sp, i=1:4, x + y <= 2*i)

macro setStageProfit(m, ex)
    m = esc(m)
    quote
        stagedata($m).stage_profit = @defExpr $(esc(gensym())) $(esc(ex))
    end
end
