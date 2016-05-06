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
- `transition`: Transition probabilities for markov chain. Either a square matrix of size `markov_states` or a vector of such matrices with length `stages`. Defaults to uniform transition probability.
- `initial_markov_state`: index of the initial markov state. If not given assumed to transition uniformly at beginning of first stage.
- `cuts_filename`: if specified, the cuts will be written to the file.
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


### Solve

#### Load previously generated cuts
If `cuts_filename` was specified in the model definition, you can load cuts via `loadcuts!(m)`, otherwise you can load cuts using `loadcuts!(m::SDDPModel, filename::ASCIIString)`.

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
  - The heuristic is run each time `frequency` cuts have been added to the model.
- `Deterministic(frequency::Int)`: solves an LP for each constraint to determine whether cut is currently defining at any point in state space.
  - The heuristic is run each time `frequency` cuts have been added to the model.

Note: Tuning the `frequency` parameter is a trade off between the model creation time, and the model solution time. If the model takes a short time to solve relative to the creation time, a low value for `frequency` may hurt performance. However, if the model takes a long time to solve relative to the creation time, aggressive cut selection (i.e.  `frequency` is small) may help performance.

#### Parallelisation
You can change the parallelisation options by constructing a `Parallel()` type and passing it one or both of the following in any order. The default is to run the pass in serial mode on a single processor.
- `ForwardPass()`
  - Forward (simulation) passes are conducted independently on all available processors. This option scales well with the number of processors due to the independent nature of estimating the objective of the policy.
- `BackwardPass(n::Int)`
 - `n` backwards (cutting) passes are computed independently on all available processors before the cuts are shared and syncronised across all processors.
 - Tuning this parameter is a trade off between the parallel overhead (copying the model and cuts between processors) and the extra cuts discovered by solving in parallel. For small models, it is likely that low values of `n` may reduce performance due to this additional overhead. In addition, since cuts are discovered independently, cuts generated on one processor will not benefit from the cuts generated on other processors until they are combined.

#### The actual solve
The `solve(m::SDDPModel [; kwargs...])` function solves the SDDP model `m`. There are the following keyword parameters:
- `maximum_iterations::Int` default = `1`
  - The maximum number of iterations (cutting passes, convergence testing) to complete before termination
- `convergence::Convergence([;simulations=1, frequency=1, terminate=false, quantile=0.95])`
  - `simulations::Int` default = `1`
    - The number of realisations to conduct when testing for convergence
  - `frequency::Int` default = `1`
    - Simulate the expected cost of the policy (using n=`simulations`) every `frequency` iterations and output to user
  - `terminate::Bool` default = `false`
    - If a convergence test is conducted with the bounds found to have converged, and `terminate=true`, method will terminate. Otherwise the method will terminate at `maximum_iterations`.
      - We choose to default this to false since if there is high variance in objective, the method may terminate earlier than desired.
  - `quantile::Float64` default = `0.95`
    - Level of confidence interval to construct when testing for convergence
- `risk_measure::RiskMeasure` default = `Expectation()`
  - See section Risk Measures above.
- `cut_selection::CutSelectionMethod` default = `NoSelection()`
  - See section Cut Selection above.
- `parallel::Parallel` default = `Parallel()`
  - See section Parallelisation above.

### Simulate Policy
The `simulate(m::SDDPModel, n::Int, variables::Vector{Symbol}; parallel::Bool=false)` function simulates `n` realisations of the policy given by a converged SDDP model `m`. It returns a dictionary with an entry for each variable given in `variables`. Each dictionary entry is a vector corresponding to the stages in the model. Each item in the vector is a vector of the `n` values that the variable took in the `n` realisations. If `parallel=true` the the `n` realisations are distributed across all available processors and computed independently.

```julia
results = simulate(m, 5, [:x])

# The five realisations of the x variable in stage 1
results[:x][1] --> [1, 1.5, 1, 1.25, 1.25]
```

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
