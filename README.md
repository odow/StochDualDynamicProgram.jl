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
- `conf_level`: defaults to `0.95`.
- `cuts_filename`: if specified, the cuts will be written to the file.
- `solver`: MathProgBase compliant solver that returns duals from a linear program

### Describing the Stage Problems
We now need to define the stage problems.

#### State Variables
We can define a new state variable in the stage problem `sp` using the `@defStateVar` macro:

It consists of three arguments

1. the stage problem model object.
2. the value of the state variable at the END of the stage. It can consist of any valid JuMP `@defVar` syntax.
3. the value of the state variable at the BEGINNING of the stage. It must consist of a keyword argument.

```julia
@defStateVar(sp, x >= 0.5, x0=1)
```
Both `x` and `x0` are JuMP variables.

Alternatively, we can use indexing just as we would in a JuMP `@defVar` macro:
```julia
X0 = [3., 2.]
@defStateVar(sp, x[i=1:2], x0=X0[i])
```
In this case, both `x` and `x0` are JuMP dicts that can be indexed with the keys `1` and `2`.
All the indices must be specified in the second argument, but they can be referred to in the third argument. The indexing of `x0` will be identical to that of `x.`

#### Stage Profit
Define the stage profit for the stage problem `sp` using the `@setStageProfit` macro. If our stage objective is `min cᵀx+θ` then:
```julia
@setStageProfit(sp, dot(c, x))
```
This can be any valid JuMP syntax. The value/cost to go is handled automatically by the solver.

#### Scenario constraints
You can add scenarios to the stage problem with the `@addScenarioConstraint` macro. It has three arguments
1. stage problem model object
2. keyword definition of scenario set. Length of the set must be identical to the value of the `scenarios` keyword in `SDDPModel()`
3. JuMP constraint. Can be any valid JuMP constraint syntax. It must include the keyword from the second argument

```julia
@addScenarioConstraint(sp, RHS=rand(3), a*x <= b + RHS)
```
The scenario sets must be indexable by the integers `1` to `scenarios` where `scenarios` is the value specified in the `SDDPModel` definition.

If you add multiple scenario constraints the length of each scenario set must be identical. Scenario One corresponds to `RHS=Ω[1,i]` for all `i`.

```julia
Ω = rand(3,4)
for i=1:4
  @addScenarioConstraint(sp, RHS=Ω[:,i], a*x[i] <= b + RHS)
end
```

#### Load previously generated cuts
If `cuts_filename` was specified in the model definition, you can load cuts via `load_cuts!(m)`, otherwise you can load cuts using `load_cuts!(m::SDDPModel, filename::ASCIIString)`.

### Solve
The `solve(m::SDDPModel [; kwargs...])` function solves the SDDP model `m`. There are the following keyword parameters:
- `simulation_passes::Int` default = `1`
  - The number of realisations to conduct when testing for convergence
- `maximum_iterations::Int` default = `1`
  - The maximum number of iterations (cutting passes, convergence testing) to complete before termination
- `convergence_test_frequency::Int` default = `1`
  - Simulate the expected cost of the policy (using n=`simulation_passes`) every `convergence_test_frequency` iterations and output to user
  - If `convergence_test_frequency=0`, never test convergence. Terminate at `maximum_iterations`
- `beta_quantile::Float64 ∈ (0, 1]` default = `1`
  - The CVar β quantile quantile for nested risk aversion
  - `beta_quantile=1` is identical to expectation
- ` risk_lambda::Float64 ∈ [0, 1]` default = `1`
  - Convex weight between Expectation and CVar (1=Expectation, 0=CVar) for nested risk aversion
  - `risk_lambda * Expectation + (1 - risk_lambda) * CVar(β)`
- `cut_selection_frequency::Int` default = `0`
  - Number of cutting passes to conduct before running cut selection algorithm
  - If `cut_selection_frequency=0` no cut deletion is conducted.
  - Tuning this parameter is a trade off between the model creation time, and the model solution time. If the model takes a short time to solve relative to the creation time, a low value for `cut_selection_frequency` may hurt performance. However, if the model takes a long time to solve relative to the creation time, aggressive cut selection (i.e.  `cut_selection_frequency` is small) may help performance.
- `cut_selection_method  ∈ {LevelOne(), Deterministic()}` default = `LevelOne()`
  - `LevelOne()`: removes those cuts that are level one dominated (de Matos, Philpott, Finardi (2015). Improving the Peformance of Stochastic Dual Dynamic Programming. Journal of Computational and Applied Mathematics 290: 196-208). That is, it keeps cuts that create the best point approximation to the value function at at least one of the sample points visited so far by the algorithm.
  - `Deterministic()`: solves an LP for each constraint to determine whether cut is currently defining at any point in state space.
  - `LevelOne` is a heuristic method and has been found to work well. However, you may want to compare the objectives with either no cut selection, or the deterministic method.
- `cuts_per_processor::Int` default = `0`
  - If `cuts_per_processor>0`, `cuts_per_processor` cuts are computed on each available processor before being collected and combined to remove duplicates. The updated set of cuts is then passed back to all processors and a new set of `cuts_per_processor` cuts are computed. In addition, convergence test (forward simulation) passes are divided up to each available processor and simulated in parallel.
  - If `cuts_per_processor=0`, method runs in serial mode.
  - Tuning this parameter is a trade off between the parallel overhead (copying the model and cuts between processors) and the extra cuts discovered by solving in parallel. For small models, it is likely that low values of `cuts_per_processor` may reduce performance due to this additional overhead. In addition, since cuts are discovered independently, cuts generated on one processor will not benefit from the cuts generated on other processors until they are combined.
- `convergence_termination::Bool` default = `false`
  - If a convergence test is conducted with the bounds found to have converged, and `convergence_termination=true`, method will terminate.
  - If this is false, the method will terminate at `maximum_iterations`
  - We choose to default this to false since if there is high variance in objective, the method may terminate earlier than desired.

### Simulate Policy
The `simulate(m::SDDPModel, n::Int, variables::Vector{Symbol}; parallel::Bool=false)` function simulates `n` realisations of the policy given by a converged SDDP model `m`. It returns a dictionary with an entry for each variable given in `variables`. Each dictionary entry is a vector corresponding to the stages in the model. Each item in the vector is a vector of the `n` values that the variable took in the `n` realisations. If `parallel=true` the the `n` realisations are distributed across all available processors and computed independently.

```julia
results = simulate(m, 5, [:x])

# The five realisations of the x variable in stage 1
results[:x][1] --> [1, 1.5, 1, 1.25, 1.25]
```
