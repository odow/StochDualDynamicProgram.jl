# StochDualDynamicProgram

<!-- [![Build Status](https://travis-ci.org/odow/StochDualDynamicProgram.jl.svg?branch=master)](https://travis-ci.org/odow/StochDualDynamicProgram.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/t32f352w4ngxappk/branch/master?svg=true)](https://ci.appveyor.com/project/odow/stochdualdynamicprogram-jl/branch/master)
[![codecov.io](https://codecov.io/github/odow/StochDualDynamicProgram.jl/coverage.svg?branch=master)](https://codecov.io/github/odow/StochDualDynamicProgram.jl?branch=master) -->

## Installation
This package is unregistered so you will need to `Pkg.clone` it as follow:
```julia
Pkg.clone("https://github.com/odow/StochDualDynamicProgram.jl.git")
```

## Usage
See the `/examples` folder for example usage but briefly

### Initialising the model object
The first step is to initialise the SDDP model object. We do this using the following syntax:

If we have more than one markov state:
```julia
m = SDDPModel([;kwargs...]) do sp, stage, markov_state
  # Stage problem definition where `sp` is a `JuMP.Model`object,
end
```


Otherwise if we have a single markov state
```julia
m = SDDPModel([;kwargs...]) do sp, stage
  # Stage problem definition
end
```

`SDDPModel` takes the following keyword arguments. The first two are essential:
+ `stages`: the number of stages in the model
+ `value_to_go_bound`: a valid bound on the initial value/cost to go. i.e. for maximisation this may be some large positive number, for minimisation this may be some large negative number.

The following arguments are optional:
- `sense`: must be either `:Max` or `:Min`. Defaults to `:Max`.
- `markov_states`: the number of markov states. Defaults to `1`.
- `scenarios`: the number of scenarios in each markov state. Defaults to `1`.
- `scenario_probability`: support vector for the scenarios. Defaults to uniform distribution.
- `transition`: Transition probabilities for markov chain. Either a square matrix of size `markov_states` or a vector of such matrices with length `stages`. Defaults to uniform transition probability.
- `initial_markov_state`: index of the initial markov state. If not given assumed to transition uniformly at beginning of first stage.
- `solver`: MathProgBase compliant solver that returns duals from a linear program

### Describing the Stage Problems
We now need to define the stage problems.

#### State Variables
We can define a new state variable in the stage problem `sp` using the `@state` macro:

It consists of three arguments

1. the stage problem model object.
2. the value of the state variable at the END of the stage. It can consist of any valid JuMP `@variable` syntax.
3. the value of the state variable at the BEGINNING of the stage. It must consist of a keyword argument.

```julia
@state(sp, x >= 0.5, x0=1)
```
Both `x` and `x0` are JuMP variables.

Alternatively, we can use indexing just as we would in a JuMP `@variable` macro:
```julia
X0 = [3., 2.]
@state(sp, x[i=1:2], x0=X0[i])
```
In this case, both `x` and `x0` are JuMP dicts that can be indexed with the keys `1` and `2`.
All the indices must be specified in the second argument, but they can be referred to in the third argument. The indexing of `x0` will be identical to that of `x.`

#### Stage Profit
Define the stage profit for the stage problem `sp` using the `@stageprofit` macro. If our stage objective is `min cᵀx+θ` then:
```julia
@stageprofit(sp, dot(c, x))
```
This can be any valid JuMP syntax. The value/cost to go is handled automatically by the solver.

#### Scenario constraints
You can add scenarios to the stage problem with the `@scenarioconstraint` macro. It has three arguments
1. stage problem model object
2. keyword definition of scenario set. Length of the set must be identical to the value of the `scenarios` keyword in `SDDPModel()`
3. JuMP constraint. Can be any valid JuMP constraint syntax. It must include the keyword from the second argument

```julia
@scenarioconstraint(sp, RHS=rand(3), a*x <= b + RHS)

RHS = rand(3)
@scenarioconstraint(sp, i=1:3, a*x <= b + c*RHS[i])
```

If you add multiple scenario constraints the length of each scenario set must be identical. Scenario One corresponds to `RHS=Ω[1,i]` for all `i`.

```julia
Ω = rand(3,4)
for i=1:4
  @scenarioconstraint(sp, RHS=Ω[:,i], a*x[i] <= b + RHS)
end
```

There is also a plural version of the macro similar to JuMP's `@constraints`:

```julia
Ω = rand(3,2)
@scenarioconstraints(sp, i=1:3, begin
  x[i] <= Ω[i, 1]

  y[i] <= Ω[i, 2]
end)
```


### Load previously generated cuts
You can load cuts from a previous solve using `loadcuts!(m::SDDPModel, filename::ASCIIString)`.

### Solve

The `solve(m::SDDPModel [; kwargs...])` function solves the SDDP model `m`. There are the following keyword parameters:
- `maximum_iterations::Int` default = `1`
  - The maximum number of iterations to complete before termination
- `policy_estimation::PolicyEstimator` default = `MonteCarloEstimator()`
  - See section Policy Estimation below
- `bound_convergence::BoundConvergence` default = `BoundConvergence()`
  - See section Bound Convergence below
- `forward_pass::ForwardPass` default = `ForwardPass()`
  - See section Forward Pass Options below.
- `backward_pass::BackwardPass` default = `BackwardPass()`
  - See section Backward Pass Options below.
- `risk_measure::RiskMeasure` default = `Expectation()`
  - See section Risk Measures below.
- `cut_selection::CutSelectionMethod` default = `NoSelection()`
  - See section Cut Selection below.
- `parallel::Parallel` default = `Serial()`
  - See section Parallelisation below.
- `output::ASCIIString` default = `nothing`
  - Prints the log trace to the file specified by `output`.
- `cut_output_file::ASCIIString` default = `nothing`
  - Save the cuts generated to the file specified by `cut_output_file`. These can be loaded later using `loadcuts!`.

#### Policy Estimation
We can control the policy estimator with the `MonteCarloEstimator([;frequency=0, minsamples=10, maxsamples=minsamples, step=0, terminate=false, confidencelevel=0.95, antitheticvariates=false])` constructor. The options are:
- `frequency::Int`: the policy estimator is run every `frequency` iterations. If `frequency < 1` then the estimator will never run.
- `min::Int`: minium number of monte-carlo simulations to conduct when estimating bound using sequential sampling.
- `max::In`: maximum number of monte-carlo simulations to conduct when estimating bound using sequential sampling.
- `step::Int`: step size for sequential sampling
- `terminate::Bool`: algorithm terminates if confidence interval from policy estimation after `max` simulations contains the lower (if minimising) bound.
- `confidencelevel::Float64`: confidence level for the confidence interval construction.
- `antitheticvariates::Bool`: sample using antithetic variates technique for variance minimisation.

#### Bound Convergence
One termination criteria for the SDDP algorithm is to terminate after the lower (if minimising) bound fails to improve by a certain tolerance after a set number of iterations. This behaviour can be controlled with the `BoundConvergence([;after=0, tol=0.])` constructor. If `after < 1` then the algorithm will not terminate due to bound convergence. Otherwise, the algorithm will terminate after `after` iterations where the bound improves by less than `tol`.

#### Forward Pass Options
You can specify forward pass options with the constructor `ForwardPass([;scenarios=1, regularisation=NoRegularisation(), importancesampling=false])`.
- `scenarios` is the number of scenarios to sample in one iteration of the forward pass. It can either by an interger, or a range/vector of integers specifying the number of samples by iteration. If there are more iterations than elements in `scenarios`, future iterations will use the last element in the list.
  - examples: `scenarios = 1:10` or `scenarios = [1, 1, 10]`.
- `regularisation` specifies the regularisation to be used on the forward pass. Valid options are
  - `NoRegularisation()`: No regularisation is conducted.
  - `LinearRegularisation([inital=1[, decay_rate=0.95]])`: The objective is penalised by the sum of the absolute values of the difference between the state variables and their value on the previous forward pass.
  - `QuadraticRegularisation([inital=1[, decay_rate=0.95]])`: The objective is penalised by the sum of squares of the difference between the state variables and their value on the previous forward pass.
  Both the `Linear` and `Quadratic` regularisation penalties are multiplied by the coefficient `initial x decay_rate ^ k` where `k` is the iteration number.
- `importancesampling=true` samples the forward pass with uniform probability instead of those given by the tranistion matrix and scenario support vector.

#### Backward Pass Options
You can specify backward pass options with the constructor `BackwardPass([;multicut=false])`. If `multicut` is `true` then a cut will be added to every markov state at every state for every forward trajectory.

#### Risk Measures
You can choose the risk measure to be applied to the model. It currently accepts
 - `Expectation()` which optimises under expectation
 - `NestedCVar([;beta=1., lambda=1.])` which optimises a convex combination of Expecation and Nested CVar (i.e. `lambda * Expectation + (1 - lambda) * CVar(beta)`)
   - `beta::Float64 ∈ (0, 1]` default = `1`
     - `beta=1` is identical to `Expectation()`
   - `lambda::Float64 ∈ [0, 1]` default = `1`

#### Cut Selection
You can choose the cut selection method to be apply to the model. It currently accepts
- `NoSelection()`
  - No cut selection is applied
- `LevelOne(frequency::Int)`
  - This heuristic removes those cuts that are level one dominated (de Matos, Philpott, Finardi (2015). Improving the Peformance of Stochastic Dual Dynamic Programming. Journal of Computational and Applied Mathematics 290: 196-208).
  - The algorithm is run every `frequency` iterations.
- `Deterministic(frequency::Int)`: solves an LP for each constraint to determine whether cut is currently defining at any point in state space.
  - The algorithm is run every `frequency` iterations.

Note: Tuning the `frequency` parameter is a trade off between the model creation time, and the model solution time. If the model takes a short time to solve relative to the creation time, a low value for `frequency` may hurt performance. However, if the model takes a long time to solve relative to the creation time, aggressive cut selection (i.e.  `frequency` is small) may help performance.

#### Parallelisation
You can change the parallelisation options with the `Parallel([;forward=true, backward=true, montecarlo=true])` constructor. The default, `Serial()`, is to run the pass in serial mode on a single processor.
- `montecarlo::Bool = true`: Convergence test (simulation) passes are conducted independently on all available processors. This option scales well with the number of processors due to the independent nature of estimating the objective of the policy.
- `forwardpass::Bool = true`: Forward pass scenarios are computed independently in parallel (when more than one scenario is sampled per iteration).
- `backwardpass::Bool = true`: Cuts are computed independently for each trajectory sampled on the forward pass.

### Simulate Policy
The `simulate(m::SDDPModel, n::Int, variables::Vector{Symbol}; parallel::Bool=false)` function simulates `n` realisations of the policy given by a converged SDDP model `m`. It returns a dictionary with an entry for each variable given in `variables`. Each dictionary entry is a vector corresponding to the stages in the model. Each item in the vector is a vector of the `n` values that the variable took in the `n` realisations. If `parallel=true` the the `n` realisations are distributed across all available processors and computed independently.

```julia
results = simulate(m, 5, [:x])

# The five realisations of the x variable in stage 1
results[:x][1] --> [1, 1.5, 1, 1.25, 1.25]
```

In addition to the variables you specify via the `variables` option, the following keys are stored

 - `:Current`:  stage cost
 - `:Future`:   value of the approximated future cost variable
 - `:Markov`:   index (1, 2, ... , M) of the markov state
 - `:Scenario`: index (1,2,...,S) of the scenario

#### Historical Simulation
Often we may wanto to evaluate the policy give a particular historical realisation. This may contain correlation structures between our random variables that are not captured by the markov chain or independent scenarios. We can do this with a few modifications to the model.

First, we need to name our `@scenarioconstraint`s. i.e.
```julia
@scenarioconstraint(sp, con_name_1, i=rand(3), x<=i)

@scenarioconstraints(sp, i=rand(3), begin
  con_name_1, x<=i
  con_name_2, i - 0.5 <= x
end)
```

Next, when calling the `simulate` function, we need to supply it the RHS values for each of the named constraints for each stage as keyword arguments. *This is different to the construction of the constraints where we could used indexed sets.* We also drop the number of simulation passes since only one historical realisation is allowed.

```julia
results = simulate(m, [:x], con_name_1=rand(10), con_name_2=rand(10))

# The one realisation of the x variable in stage 1
results[:x][1] --> [1.25]
```

#### Visualisation
It is possible to create an interactive visualisation of the simulated policy with the `@visualise` macro.

```julia
@visualise(results, (stage, replication), begin
	results[:Current][stage][replication], "Accumulated Profit (\$)", (cumulative=true)

	results[:Current][stage][replication], "Week Profit (\$)"

	results[:reservoir][stage][replication][:upper], "Upper Reservoir"
	results[:reservoir][stage][replication][:lower], "Lower Reservoir"

	Price[stage, results[:Markov][stage][replication]], "Price"
end)
```
