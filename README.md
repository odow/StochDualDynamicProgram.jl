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
 - `simulation_passes`: The number of forward simulation passes to conduct when estimating the objective
 - `maximum_iterations`: the number of cutting passes to make before testing for convergence
 - `log_frequency`: the number of cutting passes to make before testing for convergence
 - `beta_quantile`: The β quantile for CVar
 - `risk_lambda`: Convex weight between Expectation and CVar (1=Expectation, 0=CVar)

### Simulate Policy
The `simulate(m::SDDPModel, n::Int, variables::Vector{Symbol})` function simulates `n` realisations of the policy given by a converged SDDP model `m`. It returns a dictionary with an entry for each variable given in `variables`. Each dictionary entry is a vector corresponding to the stages in the model. Each item in the vector is a vector of the `n` values that the variable took in the `n` realisations.

```julia
results = simulate(m, 5, [:x])

# The five realisations of the x variable in stage 1
results[:x][1] --> [1, 1.5, 1, 1.25, 1.25]
```
