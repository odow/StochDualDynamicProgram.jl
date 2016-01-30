# StochDualDynamicProgram

[![Build Status](https://travis-ci.org/odow/StochDualDynamicProgram.jl.svg?branch=master)](https://travis-ci.org/odow/StochDualDynamicProgram.jl)

## Installation
This package is unregistered so you will need to `Pkg.clone` it as follow:
```julia
Pkg.clone("https://github.com/odow/StochDualDynamicProgram.jl.git")
```

## Usage
See the `/examples/hydro.jl` file for example usage. But briefly:

### `SDDPModel()`

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
