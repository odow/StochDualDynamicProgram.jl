function State(m::Model, lower::Number, upper::Number, name)
    Variable(m,lower,upper,:Cont,utf8(string(name)),NaN)
end

function State0(m::Model, init::Number, name, name0)
    v0 = Variable(m,-Inf,Inf,:Cont,utf8(string(name0)),init)
    c = @addConstraint(m, v0 == init)
    if haskey(m.ext[:dual_constraints], name)
        push!(m.ext[:dual_constraints][name], c)
    else
        m.ext[:dual_constraints][name] = Any[c]
    end
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
    x0_value = x0.args[2]

    quotx0name = quot(getname(x0_name))
    escx0name = esc(getname(x0_name))

    quotvarname = quot(getname(var))
    escvarname  = esc(getname(var))

    if isa(var,Symbol)
        # Easy case - a single variable
        return assert_validmodel(m, quote
            $(esc(var)) = State($m, $lb, $ub, $quotvarname)
            registervar($m, $quotvarname, $escvarname)
            push!($m.ext[:state_vars], $quotvarname)

            $(esc(x0_name)) = State0($m, $x0_value, $quotvarname, $quotx0name)
            registervar($m, $quotx0name, $escx0name)

            $m.ext[:dual_constraints][$quotvarname] = [@addConstraint(m, $(esc(x0_name)) == $x0_value)]

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
        push!($m.ext[:state_vars], $quotvarname)
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

function StageProblem(scenarios::Int=1)
    sp = Model()
    sp.ext[:is_sp] = true
    sp.ext[:state_vars] = Any[]
    sp.ext[:dual_constraints] = Dict{Any, Vector{ConstraintRef}}()
    sp.ext[:theta] = nothing
    sp.ext[:old_scenario] = (0,0)
    sp.ext[:scenario_constraints] = Tuple{Any, Vector{Any}}[]
    sp.ext[:LastScenario] = 0
    sp.ext[:CurrentScenario] = 0
    sp.ext[:objective_value] = zeros(scenarios)
    sp.ext[:dual_values] = Dict{Any, Vector{Float64}}()
    return sp
end

sp = StageProblem()
@defStateVar(sp, 0 <= x <= 1, y0=0.5)

sp.ext
m=Model();m.ext[:state_vars] = Symbol[];m.ext[:dual_constraints] = Dict{Any, Any}()
@JuMPdefStateVar(m, x>=0, x0=1.12)
solve(m)
getValue(x0)
m.ext



m=Model();m.ext[:state_vars] = Symbol[];m.ext[:dual_constraints] = Dict{Symbol, Any}()

S = [:upper, :lower]
A = Dict{Symbol, Float64}(:upper=>2., :lower=>1.)

@defStateVar(m, x[r=S]>=A[r], x0=A[r])
m.ext
m.colLower
solve(m)
getValue(x0[:lower])


println("The End.")
