type StageData
    state_vars::Vector{JuMP.Variable}
    dual_constraints::Vector{JuMP.ConstraintRef}
    theta::Union{Void, JuMP.Variable}
    old_scenario::Tuple{Int, Int}
    scenario_constraints::Vector{Tuple{Any, Vector{Any}}}
    scenario_constraint_names::Dict{Symbol, Int}
    last_scenario::Int
    current_scenario::Int
    objective_value::Vector{Float64}
    dual_values::Vector{Vector{Float64}}
    stage_profit
    cut_data::Dict{Tuple{Float64, Float64}, Vector{Vector{Float64}}}
end
StageData(scenarios::Int=1) = StageData(Variable[], ConstraintRef[], nothing, (0,0), Tuple{Any, Vector{Any}}[], Dict{Symbol, Int}(), 0, 0, zeros(scenarios), Vector{Float64}[], nothing, Dict{Tuple{Float64, Float64}, Vector{Vector{Float64}}}())

"""
Instaniates a new StageProblem which is a JuMP.Model object with an extension
dictionary.
"""
function StageProblem(scenarios::Int=1)
    sp = Model()
    sp.ext[:data] = StageData(scenarios)
    return sp
end

function stagedata(m::Model)
    @assert haskey(m.ext, :data)
    return m.ext[:data]::StageData
end

function Base.copy(sd::StageData)
    return StageData(sd.state_vars, sd.dual_constraints, sd.theta, sd.old_scenario, sd.scenario_constraints, sd.scenario_constraint_names, sd.last_scenario, sd.current_scenario, sd.objective_value, sd.dual_values, sd.stage_profit, sd.cut_data)
end

# Check this is a StageProblem
is_sp(m::JuMP.Model) = isa(stagedata(m), StageData)
is_sp(m) = false

type SDDPModel{T<:Union{Array{Float64,2}, Vector{Array{Float64, 2}}}, M<:Real}
    stages::Int
    markov_states::Int
    scenarios::Int
    sense::Union{Type{Val{:Min}}, Type{Val{:Max}}}

    stage_problems::Array{JuMP.Model}
    transition::T#Array{Union{Float64, Array{Float64, 2}}, T}
    scenario_probability::WeightVec
    init_markov_state::Int
    confidence_interval::Tuple{Float64, Float64}
    valid_bound::Float64
    QUANTILE::Float64
    LPSOLVER::MathProgBase.AbstractMathProgSolver
    value_to_go_bound::M
    beta_quantile::Float64
    risk_lambda::Float64
    weightings_matrix::Array{Float64, 2}
    cuts_filename::Union{Void, ASCIIString}

    build_function!::Function # function that builds the stage problems
    stagecuts::Array{StageCuts, 2}
end

"""
Creates an SDDPModel type

    SDDPModel(kwargs...)

Keyword Arguments:
cuts_filename     - ASCIIString filename for cut output. If specified cuts will be written.
initial_markov_state - the scenario at time 0. Model transitions at the start of time period 1
markov_states     - the number of markov states
scenarios         - the number of stagewise independent scenarios in a markov state
scenario_probability - vector with probability of each scenario occuring
sense             - :Max or :Min
solver            - AbstractMathProgBase solver capable of returning dual variables. Defaults to Clp.
stages            - the number of stages
transition        - the transition matrix. Can be N*N where N is the number of scenarios, or M*N*N where M is the number of stages.
value_to_go_bound - Initial bound on value to go
"""
function SDDPModel(;
    cuts_filename=nothing,
    initial_markov_state=0,
    markov_states=1,
    scenarios=1,
    scenario_probability=nothing,
    sense=:Max,
    solver=ClpSolver(),
    stages=1,
    transition=nothing,
    value_to_go_bound=1000.,
    build_function = ()->()
    )

    # Check non-zero stages, markov states and scenarios
    @assert stages >= 1
    @assert markov_states >= 1
    @assert scenarios >= 1

    if !(sense in [:Max, :Min])
        error("Sense must be either :Max or :Min. Currently = $(sense).")
    end

    if transition==nothing
        # User hasn't specified transition matrix. Assume uniform
        transition = 1 / markov_states * ones(markov_states, markov_states)
    end
    T = length(size(transition))
    if T==2
        # Stage independent transition probabilities
        @assert is_valid_transition_matrix(transition, markov_states)
    elseif T==1
        # Vector of transition matrices
        @assert length(transition) == stages
        for i=1:stages
            @assert is_valid_transition_matrix(transition[i], markov_states)
        end
    end
    # TODO :: Case where Transition matrix is Array{Float64, 3}

    if scenario_probability == nothing
        scenario_probability = ones(scenarios) / scenarios
    else
        @assert abs(sum(scenario_probability) - 1) < 1e-5
    end

    my_inf = (sense==:Max?Inf:-Inf)
    sense_type = (sense==:Max?Val{:Max}:Val{:Min})

    SDDPModel(stages,markov_states,scenarios, sense_type, Array(JuMP.Model, (stages, markov_states)), transition, WeightVec(scenario_probability), initial_markov_state, (-my_inf, -my_inf), my_inf, 0.95, solver, value_to_go_bound, 1., 1., zeros(markov_states, scenarios), cuts_filename, build_function, Array(StageCuts, (0, 0)))
end

"""
This fuction checks that the matrix T is a valid markov transition matrix
"""
function is_valid_transition_matrix(T::Array{Float64, 2}, n::Int)
    if size(T)[1] != size(T)[2]
        error("Transition matrix must be square")
    end
    @assert size(T)[1] == n
    for i=1:n
        if abs(sum(T[i,:]) - 1.) > 1e-10
            return false
        end
    end
    return true
end

# This will probably break at some point
# https://groups.google.com/forum/#!searchin/julia-users/arguments$20anonymous$20function/julia-users/QcgdNZd-sI8/dgIYSxozZzUJ
if VERSION > v"0.5-"
    arglength(f::Function)=length(Base.uncompressed_ast(methods(f).defs.func).args[1])-1
else
    arglength(f::Function)=length(Base.uncompressed_ast(f.code.def).args[1])
end

function is_zero_objective{T}(ex::JuMP.GenericQuadExpr{T, JuMP.Variable})
    return length(ex.qcoeffs) == 0 && length(ex.qvars1) == 0 && length(ex.qvars2) == 0 && is_zero_objective(ex.aff)
end
function is_zero_objective{T}(ex::JuMP.GenericAffExpr{T, JuMP.Variable})
    return length(ex.coeffs) == 0 && length(ex.vars) == 0 && ex.constant == 0.0
end

function SDDPModel(
    build_subproblem!::Function;
    cuts_filename=nothing,
    initial_markov_state=0,
    markov_states=1,
    scenarios=1,
    scenario_probability=nothing,
    sense=:Max,
    stages=1,
    transition=nothing,
    solver=ClpSolver(),
    value_to_go_bound=NaN
    )
    if isnan(value_to_go_bound)
        error("You must specify the option [value_to_go_bound] when creating an SDDPModel.")
    end

    # Initialise SDDPModel object
    m = SDDPModel(
        sense=sense,
        stages=stages,
        markov_states=markov_states,
        scenarios=scenarios,
        scenario_probability=scenario_probability,
        transition=transition,
        initial_markov_state=initial_markov_state,
        solver=solver,
        value_to_go_bound=value_to_go_bound,
        cuts_filename=cuts_filename,
        build_function=build_subproblem!
    )
    create_subproblems!(m)

    return m
end

function create_subproblems!(m::SDDPModel)
    # For every stage
    for stage=1:m.stages
        # For every markov state
        for markov_state=1:m.markov_states
            # Create a new stage problem with [scenarios] number of scenarios
            sp = StageProblem(m.scenarios)

            # Set the solver
            setsolver(sp, m.LPSOLVER)

            # User has specified markov states
            if arglength(m.build_function!)==3
                m.build_function!(sp, stage, markov_state)

            # No markov states specified
            elseif arglength(m.build_function!)==2
                # Double check they only mean one markov state
                @assert m.markov_states == 1

                m.build_function!(sp, stage)
            else
                error("Invalid number of arguments")
            end

            # Initialise storage for scenario duals now we know the number of state variables
            for i in 1:length(stagedata(sp).state_vars)
                push!(stagedata(sp).dual_values, zeros(m.scenarios))
            end

            # If the user hasn't specified an objective
            if is_zero_objective(getobjective(sp))
                if stage==m.stages
                    # If its the last stage then its just the stage profit
                    set_objective!(m.sense, sp)
                else
                    # Otherwise create a value/cost to go variable
                    set_objective!(m.sense, sp, m.value_to_go_bound)
                end
            end
            # Store the stage problem
            m.stage_problems[stage, markov_state] = sp

        end
    end

    return
end

function set_objective!(::Type{Val{:Min}}, sp::Model, value_to_go_bound)
    stagedata(sp).theta = @variable(sp, theta >= value_to_go_bound)
    @objective(sp, Min, stagedata(sp).stage_profit + stagedata(sp).theta)
end
function set_objective!(::Type{Val{:Max}}, sp::Model, value_to_go_bound)
    stagedata(sp).theta = @variable(sp, theta <= value_to_go_bound)
    @objective(sp, Max, stagedata(sp).stage_profit + stagedata(sp).theta)
end
function set_objective!(::Type{Val{:Min}}, sp::Model)
    @objective(sp, Min, stagedata(sp).stage_profit)
end
function set_objective!(::Type{Val{:Max}}, sp::Model)
    @objective(sp, Max, stagedata(sp).stage_profit)
end

"""
So we can copy an SDDPModel
"""
function Base.copy(m::SDDPModel)
    SDDPModel(m.stages, m.markov_states, m.scenarios, m.sense, deepcopy(m.stage_problems), copy(m.transition), m.scenario_probability, m.init_markov_state, m.confidence_interval, m.valid_bound, m.QUANTILE, m.LPSOLVER, m.value_to_go_bound, m.beta_quantile, m.risk_lambda, m.weightings_matrix, m.cuts_filename, m.build_function!, m.stagecuts)
end

"""
Transition markov states

Inputs:
m            - The SDDPModel object
stage        - the current stage
markov_state - the current markov state

Returns:
    the new markov state::Int
"""
function transition(m::SDDPModel, stage::Int, markov_state::Int, r::Float64)
    for new_markov_state=1:m.markov_states
        k = get_transition(m, stage, markov_state, new_markov_state)
        if r <= k
            return new_markov_state
        else
            r -= k
        end
    end
    return m.markov_states
end
transition(m::SDDPModel, stage::Int, markov_state::Int) = transition(m, stage, markov_state, rand())

function get_transition(m::SDDPModel{Array{Float64, 2}}, stage::Int, current_markov_state::Int, new_markov_state::Int)
    if current_markov_state==0
        return 1./m.markov_states
    else
        return m.transition[current_markov_state,new_markov_state]::Float64
    end
end
function get_transition(m::SDDPModel{Array{Array{Float64, 2}, 1}}, stage::Int, current_markov_state::Int, new_markov_state::Int)
    if current_markov_state==0
        1./m.markov_states
    else
        m.transition[stage][current_markov_state,new_markov_state]::Float64
    end
end

"""
This function loads cuts from a file.
"""
function loadcuts!(m::SDDPModel, filename::ASCIIString)
    open(filename, "r") do f
        while true
            line = readline(f)
            (line == nothing || line == "") && break
            line = split(strip(line), ",")
            @assert length(line) >= 3
            stage = parse(Int, line[1])
            markov_state = parse(Int, line[2])
            sp = m.stage_problems[stage, markov_state]
            theta = parse(Float64, line[3])

            @assert length(line) == (3 + length(stagedata(sp).state_vars))
            if length(line) > 4
                xcoeff = map(x->parse(Float64, x), line[4:end])
            else
                xcoeff = [parse(Float64, line[4])]
            end
            @constraint(sp, stagedata(sp).theta <= theta + sum{xcoeff[i]*getvariable(sp, v), (i, v) in enumerate(stagedata(sp).state_vars)})
        end
    end
end
function loadcuts!(m::SDDPModel)
    if m.cuts_filename == nothing
        error("Please specify a file to load cuts from.")
    end
    loadcuts!(m, m.cuts_filename)
end
