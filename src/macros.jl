# ------------------------------------------------------------------------------
#
# This file contains JuMP extension macros
#
# Much of the functionality in @state is shamelessly copied from
# both the JuMP @variable

#  Copyright 2016, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import JuMP: assert_validmodel, validmodel, esc_nonconstant
import JuMP: getloopedcode, buildrefsets, getname, registervar
import JuMP: storecontainerdata, isdependent, JuMPContainerData, pushmeta!, JuMPContainer
using Base.Meta

function State(m::Model, lower::Number, upper::Number, name)
    v = Variable(m,lower,upper,:Cont,utf8(string(name)),0.5*(lower+upper))
    push!(stagedata(m).state_vars, v)
    return v
end

function State0(m::Model, init, name, name0)
    v0 = Variable(m,-Inf,Inf,:Cont,utf8(string(name0)),init)
    push!(stagedata(m).dual_constraints, @constraint(m, v0 == init))
    return v0
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
macro state(args...)
    #
    # macro state(sp, x, x0)
    #     code = quote
    #         push!(stagedata(sp).state_vars, macrocall(:variable, sp, x))
    #     end
    #     buildrefsets etc. ala https://github.com/mlubin/JuMPChance.jl/blob/master/src/macros.jl
    # end
    # 
    length(args) <= 1 &&
        error("in @state: expected model as first argument, then variable information.")
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
                    error("in @state ($var): expected <= operator after variable name.")
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
            error("in @state ($(string(x))): use the form lb <= ... <= ub.")
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
    isa(var,Expr) || error("in @state: expected $var to be a variable name")

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
