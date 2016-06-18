nonzero(x) = abs(x) > 1e-10

function fixedmodel(m::Model)
    !m.internalModelLoaded && error("Unable to fix an unsolved model")
    f = copy(m)
    # setsolver(f, ClpSolver())
    # Follow what Gurobi does for Special Ordered Sets
    for sos in m.sosconstr
        cols = Int[v.col for v in sos.terms]
        if sos.sostype == :SOS1
            for c in cols
                if nonzero(m.colVal[c])
                    @constraint(f, sum{Variable(f, i), i=cols; i!=c} == 0)
                    break
                end
            end
        elseif sos.sostype == :SOS2
            colorder = sortperm(sos.weights)
            for (i, sosidx) in enumerate(colorder)
                c = cols[sosidx]
                if nonzero(m.colVal[c]) && i != length(colorder)
                    @constraint(f, sum{Variable(f, j), j=cols; j!=c && j != cols[colorder[i+1]]} == 0)
                    break
                elseif nonzero(m.colVal[c]) && i == length(colorder)
                    @constraint(f, sum{Variable(f, j), j=cols; j!=c && j != cols[colorder[i-1]]} == 0)
                    break
                end
            end
        else
            error("Unable to create fixed model with SOS of type $(sos.sostype)")
        end
    end
    f.sosconstr = JuMP.SOSConstraint[]

    # Fix all binary and integer variables to their current value
    for (i, v) in enumerate(m.colCat)
        if v == :Bin || v == :Int
            setcategory(Variable(f, i), :Cont)
            @constraint(f, Variable(f, i) == m.colVal[i])
        end
    end
    f
end

function solvenoncontinuous!(sp::Model, scenario::Int)
    f = fixedmodel(sp)
    status = solve(f)
    @assert status == :Optimal
    sp.linconstrDuals = f.linconstrDuals
    # store the objective value
    stagedata(sp).objective_values[scenario] = getobjectivevalue(f)
    # store the dual value for each of the state variables
    for i in 1:length(stagedata(sp).state_vars)
        stagedata(sp).dual_values[scenario][i] = f.linconstrDuals[stagedata(sp).dual_constraints[i].idx]
    end
end

function SOSII!(m::JuMP.Model, f::Function, x, lb::Float64, ub::Float64, n::Int64=21)
    #   This function adds a SOS of Type II to the JuMP model
    #      Σλ = 1         # Convexity
    #      a∙λ = z        # x value
    #      f(z) = f(a)∙λ  # approximated function value
    #      where λ is a SOS of Type II
    xx = linspace(lb, ub, n) # Set x value
    fx = f(xx)               # Estimate function value
    SOSII!(m, x, xx, fx)
end

function SOSII!(m::JuMP.Model, x, xx::Vector, fx::Vector)
    n = length(xx)
    @assert length(fx) == n

    @variables(m, begin
        y                # Estimated function value variable
        0 <= λ[1:n] <= 1 # SOS variable set
    end)
    @constraints(m, begin
        sum(λ) == 1          # Convexity
        x      == dot(xx, λ)
        y      == dot(fx, λ)
    end)
    addSOS2(m, [i*λ[i] for i =1:n])
    y  # Return estimated function value variable
end

SOSII!{T<:Real}(m::JuMP.Model, x, xx::Vector{NTuple{2, T}}) = SOSII!(m, x, [y[1] for y in xx], [y[2] for y in xx])

function bilinear(m::JuMP.Model, x::Variable, y::Variable)
    # This macro is used to replace a bilinear term in a JuMP constraint with the expansion
    #   xy = 0.25 * [(x+y)^2 - (x-y)^2]
    #     where (x+y)^2 and (x-y)^2 are approximated with SOS2

    # (x + y)^2 SOS2
    v1 = SOSII!(m, (a) -> a.^2, x + y, getLower(x) + getLower(y), getUpper(x) + getUpper(y))

    # (x + y)^2 SOS2
    v2 = SOSII!(m, (a) -> a.^2, x - y, getLower(x) - getUpper(y), getUpper(x) - getLower(y))

    # Return expression to be placed in constraint
    (0.25 * (v1 - v2))
end
