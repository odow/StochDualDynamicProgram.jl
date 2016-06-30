#
#   This example demonstrates how the user can incorporate specialised solve routines
#   into the normal stage problem solve via the JuMP `solvehook` functionality.
#
using StochDualDynamicProgram, JuMP

m = SDDPModel(stages=2, scenarios=5, value_to_go_bound=100) do sp, stage

    @state(sp, x, x0=0)

    @variables(sp, begin
        noise
        slack
    end)

    @scenarioconstraint(sp, err=rand(5), noise == err - 0.5)

    @constraint(sp, x == x0 + noise + slack)

    sp.ext[:penalty_constraint] = @constraint(sp, slack == 0)

    @stageprofit(sp, x)

    function twostagesolve(sp; kwargs...)
        status = solve(sp; ignore_solve_hook=true)
        xvalue = getvalue(getvariable(sp, :x))
        if xvalue < 0
            JuMP.setRHS(sp.ext[:penalty_constraint], abs(xvalue))
        end
        status = solve(sp; ignore_solve_hook=true)
        JuMP.setRHS(sp.ext[:penalty_constraint], 0.)
        return status
    end

    JuMP.setsolvehook(sp, twostagesolve)

end

results = solve(m, maximum_iterations=20)

simulations = simulate(m, 100, [:x])

for t=1:2
    for s=1:100
        @assert simulations[:x][t][s] >= 0
    end
end
