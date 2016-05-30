function backwardpass!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, riskmeasure::RiskMeasure, regularisation::Regularisation=NoRegularisation())
    for t=(T-1):-1:1
        for pass = 1:getn(m)
            setrhs!(m, pass, t)
            solveall!(m, t+1, regularisation)
            reweightscenarios!(m, t, getmarkov(m, pass, t), riskmeasure.beta, riskmeasure.lambda)
            addcut!(m, pass, t, getmarkov(m, pass, t))
        end
    end
end

# a wrapper for backward pass timings
function backwardpass!(log::SolutionLog, m::SDDPModel, nscenarios, risk_measure, regularisation, isparallel::Bool)
    tic()
    if isparallel
        parallelbackwardpass!(m, risk_measure, regularisation)
    else
        backwardpass!(m, risk_measure, regularisation)
    end
    log.cuts += nscenarios
    log.time_backwards += toq()
    return
end

function addcut!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, pass::Int, t::Int, i::Int)
    # calculate coefficients
    xbar = getx(m, pass, t)
    cut = Cut(length(xbar))
    sd1 = stagedata(m, t, i)
    for j=1:M
        for s=1:S
            sd = stagedata(m, t+1, j)
            cut.intercept += sd1.weightings_matrix[j,s] * sd.objective_values[s]
            for v = 1:length(xbar)
                cut.intercept -= sd1.weightings_matrix[j,s] * sd.dual_values[v][s] * xbar[v]
                cut.coefficients[v] += sd1.weightings_matrix[j,s] * sd.dual_values[v][s]
            end
        end
    end
    # add to cutselection storage
    add_cut!(X, stagecut(m, t, i), cut)
    # add to problem
    addcut!(X, subproblem(m, t, i), cut)

    return cut
end

function addcut!(::Type{Max}, sp::Model, cut::Cut)
    sd = stagedata(sp)
    @constraint(sp, sd.theta <= cut.intercept + dot(cut.coefficients, sd.state_vars))
end
function addcut!(::Type{Min}, sp::Model, cut::Cut)
    sd = stagedata(sp)
    @constraint(sp, sd.theta >= cut.intercept + dot(cut.coefficients, sd.state_vars))
end

# set the rhs of the t+1 subproblems using value from forward pass
function setrhs!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, pass::Int, t::Int)
    for i=1:M
        sp = subproblem(m, t+1,i)
        for j in 1:length(stagedata(sp).dual_constraints)
            JuMP.setRHS(stagedata(sp).dual_constraints[j], getx(m, pass, t, j))
        end
    end
end

# solve all the t subproblems
function solveall!{T, M, S, X, TM}(m::SDDPModel{T, M, S, X, TM}, t::Int, regularisation::Regularisation)
    for i in 1:M
        set_nonregularised_objective!(regularisation, X, subproblem(m,t,i))
        for s=1:S
            load_scenario!(subproblem(m,t,i), s)
            backsolve!(subproblem(m,t,i), s)
        end
    end
end

# Solve the subproblem sp in the backward pass storing the objective and dual coefficients
function backsolve!(sp::Model, scenario::Int)
    @assert issubproblem(sp)
    status = solve(sp)
    # Catch case where we aren't optimal
    if status != :Optimal
        sp.internalModelLoaded = false
        status = solve(sp)
        if status != :Optimal
            error("SDDP Subproblems must be feasible. Current status: $(status). I tried rebuilding from the JuMP model but it didn't work...")
        end
    end
    # store the objective value
    stagedata(sp).objective_values[scenario] = getobjectivevalue(sp)
    # store the dual value for each of the state variables
    for i in 1:length(stagedata(sp).state_vars)
        stagedata(sp).dual_values[i][scenario] = getdual(stagedata(sp).dual_constraints[i])
    end
end
