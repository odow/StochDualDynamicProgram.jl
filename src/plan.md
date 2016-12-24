# Plan of attack

SDDPModel Type
    - Stores a matrix of subproblems
        - subproblem::JuMP.Model with extension that stores
            - list of state variables. State variable has
                - a JuMP.variable
                - a JuMP.constraint for the dummy dual
            - a weighted vector of RHS constraints
                - each scenario has a list of constraint/weight pairs
            - price scenario information
                - position of each rib (abstract vector)
                - noises (abstract vector)
                - dynamics (function)
                - objective function
                - value
            - de Matos cut selection information
                - list of states visited
                - for each rib
                    - list of cut coefficients
                    - number of points where cut is dominant
                    - list of indices for each recorded state of cut which is dominant at that point
    -       - a risk measure
                - dispatches on constructcut(riskmeasure::AbstractRiskMeasure, objectives::Vector{Float64}, values::Vector{Vector{Float64}}, duals::Vector{Vector{Float64}})
                - riskmeasure is the instance of the type so can contain lots of information
                - examples include
                    - ::Expectation
                    - ::NestedCVaR
                    - ::Robust

macros
    - @state(s)
    - @scenario(s)
    - @stageobjective

    - pricestate!(sp,              # subproblem
        discretisation = 0:10,     # outgoing price
        initial        = 0.5,      # incoming price
        noise          = rand(10), # noises
        (p, w) -> (p + w),         # dynamics
        (p) -> (p * x)             # Objective can use p0 as a parameter, x as a variable
    )

solve options
    - scenario incrementation
    - max iterations
    - bound stalling
    - confidence interval
    - uniform sampling
    - expected value iterations
    - non-convexity
    - add cuts to all markov states backward pass
    - risk measure
        - see AbstractRiskMeasure above
    - parallelism
        - asyncronous cutting passes
        - simulation

what not to include
    - regularisation
