# StochDualDynamicProgram

[![Build Status](https://travis-ci.org/odow/StochDualDynamicProgram.jl.svg?branch=master)](https://travis-ci.org/odow/StochDualDynamicProgram.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/t32f352w4ngxappk/branch/master?svg=true)](https://ci.appveyor.com/project/odow/stochdualdynamicprogram-jl/branch/master)
[![codecov.io](https://codecov.io/github/odow/StochDualDynamicProgram.jl/coverage.svg?branch=master)](https://codecov.io/github/odow/StochDualDynamicProgram.jl?branch=master)

## Installation
This package is unregistered so you will need to `Pkg.clone` it as follow:
```julia
Pkg.clone("https://github.com/odow/StochDualDynamicProgram.jl.git")
```

## Usage
See the `/examples/hydro.jl` file for example usage. But briefly:

### `SDDPModel`

This function creates a new SDDP model object. `SDDPModel()` takes the following keyword arguments:

| Keyword         | Description                                      | Default     |
| ----------------| ------------------------------------------------ | ----------- |
| `sense`         | `:Max` if maximisation or `:Min` if minimisation | `:Max`      |
| `stages`        | Number of stages in the problem                  | 1           |
| `markov_states` | Number of markov states in the problem           | 1           |
| `transition`    | Transition matrix for markov states. Must be either a square matrix of size `markov_states` (same transition probabilities for each stage) or a vector of length `stages` containing square matrices of size `markov_states` (different transition probabilities in each stage). If no transition matrix is specified, then it is assumed there is an equiprobable transition matrix.  | nothing     |
| `initial_markov_state` | Initial markov state. Set to 0 if there is an equal probability of arriving at an initial state | 0 |
| `conf_level`    |  Confidence level to use when determining if bound is within simulated objective | 0.95        |
| `lpsolver` |  MathProgBase solver to use. Must return dual variables.          | ClpSolver() |
| `abs_tol`  | Absolute tolerance to use when checking if the bound is within simulated confidence interval. | `1e-8`      |
| `rel_tol`  | Relative tolerance to use when checking if the bound is within simulated confidence interval. | `1e-8`      |

#### Examples
```julia
m = SDDPModel(stages=3, markov_states=2, transition=[0.2 0.8; 0.7 0.3])
m = SDDPModel(sense=:Min, stages=3, markov_states=2, conf_level=0.99)
```
### `addStageProblem!`
This function defines a new stage problem in the SDDP model `m`.
```julia
addStageProblem!(m::SDDPModel, stage::Int, markov_state::Int) do sp
  # Stage problem definition
end
```

### `@defStateVar`
Define a new state variable in the stage problem `sp`.

#### Example
```julia
@defStateVar(sp, x, x0==1)
@defStateVar(sp, x >= 0.5, x0==1)
@defStateVar(sp, 0 <= x <= 0.5, x0==1)
```

### `@defValueToGo`
Define the value to go variable for the stage problem `sp`.

The value to go variable must be defined with some large, non-infinite bound.

If the SDDP model has the sense `:Max`, it is expected that theta optimises towards `+inf`. 
If the SDDP model has the sense `:Min`, it is expected that theta optimises towards `-inf`. 

#### Example
```julia
@defValueToGo(sp, ϑ <= 1000)
```

### `solve(m::SDDPModel [; forward_passes=1, backward_passes=1)`
This function solves the SDDP model `m`. There are two optional keyword parameters:
 - `forward_passes`: The number of forward simulation passes to conduct when estimating the objective
 - `backward_passes`: the number of cutting passes to make before testing for convergence

### `simulate(m::SDDPModel, n::Int, variables::Vector{Symbol})`
This function simulates `n` realisations of the policy given by a converged SDDP model `m`. It returns a dictionary with an entry for each variable given in `variables`. Each dictionary entry is a vector corresponding to the stages in the model. Each item in the vector is a vector of the `n` values that the variable took in the `n` realisations.

```julia
results = simulate(m, 5, [:x])

# The five realisations of the x variable in stage 1 
results[:x][1] --> [1, 1.5, 1, 1.25, 1.25]
```
