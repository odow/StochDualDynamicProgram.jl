nonzero(x) = abs(x) > 1e-10

function fixvariable!(m::Model, idx, value)
    setupperbound(Variable(m, idx), value)
    setlowerbound(Variable(m, idx), value)
end

function fixmodel!(m::Model)
    !m.internalModelLoaded && error("Unable to fix an unsolved model")
    # Follow what Gurobi does for Special Ordered Sets
    if length(m.sosconstr) > 0
        for sos in m.sosconstr
            cols = Int[v.col for v in sos.terms]
            if sos.sostype == :SOS1
                for c in cols
                    nonzero(m.colVal[c]) || fixvariable!(m, c, 0.)
                end
            elseif sos.sostype == :SOS2
                colorder = sortperm(sos.weights)
                for (i, sosidx) in enumerate(colorder)
                    c = cols[sosidx]
                    nonzero(m.colVal[c]) || fixvariable!(m, c, 0.)
                end
            else
                error("Unable to create fixed model with SOS of type $(sos.sostype)")
            end
        end
        m.sosconstr = JuMP.SOSConstraint[]
        # Do this for now since no efficient way to update/remove SOS constraints
        m.internalModelLoaded = false
    end
    # Fix all binary and integer variables to their current value
    for (i, v) in enumerate(m.colCat)
        if v == :Bin || v == :Int
            setcategory(Variable(m, i), :Cont)
            fixvariable!(m, i, m.colVal[i])
        end
    end
end

function solvenoncontinuous!(sp::Model, scenario::Int)
    lb = copy(sp.colLower)
    ub = copy(sp.colUpper)
    cat = copy(sp.colCat)
    sos = copy(sp.sosconstr)

    fixmodel!(sp)
    status = solve(sp)
    @assert status == :Optimal

    storecontinuous!(sp, scenario)

    sp.colLower = lb
    sp.colUpper = ub
    sp.colCat = cat
    sp.sosconstr = sos

end

function SOSII!(m::JuMP.Model, f::Function, x, lb::Float64, ub::Float64, n::Int64=10)
    #   This function adds a SOS of Type II to the JuMP model
    #      Σλ = 1         # Convexity
    #      a∙λ = z        # x value
    #      f(z) = f(a)∙λ  # approximated function value
    #      where λ is a SOS of Type II
    xx = linspace(lb, ub, n) # Set x value
    fx = f(xx)               # Estimate function value
    SOSII!(m, x, xx, fx)
end

function SOSII!(m::JuMP.Model, x, xx::AbstractVector, fx::AbstractVector)
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

function bilinear(x::Variable, y::Variable, n::Int=10)
    # This macro is used to replace a bilinear term in a JuMP constraint with the expansion
    #   xy = 0.25 * [(x+y)^2 - (x-y)^2]
    #     where (x+y)^2 and (x-y)^2 are approximated with SOS2
    m = x.m

    # (x + y)^2 SOS2
    v1 = SOSII!(m, (a) -> a.^2, x + y, getLower(x) + getLower(y), getUpper(x) + getUpper(y), n)

    # (x + y)^2 SOS2
    v2 = SOSII!(m, (a) -> a.^2, x - y, getLower(x) - getUpper(y), getUpper(x) - getLower(y), n)

    # Return expression to be placed in constraint
    (0.25 * (v1 - v2))
end
