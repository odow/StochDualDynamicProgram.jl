#  Copyright 2016, Oscar Dowson

using JuMP, StochDualDynamicProgram

function AR1_Process(α, ϵ, initial, N)
    inflow_min = 1e-2
    inflow_max = 1000.
    xpoints = logspace(log(inflow_min), log(inflow_max), N)
    ypoints = log(xpoints)

    m = Model()

    @variables(m, begin
        inflow_initial
        inflow

        log_inflow_initial
        log_inflow
    end)

    @constraints(m, begin
        # Initial value
        inflow_initial      == initial

        # Set up log(x) constraints
        log_inflow_initial == SOSII!(m, inflow_initial, xpoints, ypoints)
        log_inflow         == SOSII!(m, inflow,         xpoints, ypoints)

        # Dynamics constraint
        log_inflow         == α * log_inflow_initial + log(ϵ)
    end)

    solve(m)
    println("Initial Inflow")
    println("Actual:       $(getvalue(inflow_initial))")
    println("Log Estimate: $(getvalue(log_inflow_initial))")
    println("Log Actual:   $(log(getvalue(inflow_initial)))")
    println()
    println("Inflow")
    println("Actual:       $(getvalue(inflow))")
    println("Log Estimate: $(getvalue(log_inflow))")
    println("Log Actual:   $(log(getvalue(inflow)))")
    println()
    println("Inflow = ϵ * (Initial Inflow)ᵅ")
    println("Actual:    $(ϵ * getvalue(inflow_initial) ^ α)")
    println("Estimated: $(getvalue(inflow))")

end

#
# AR1_Process(α, ϵ, initial inflow, number breakpoints)
#
AR1_Process(0.5, 1.01, 200, 50)
