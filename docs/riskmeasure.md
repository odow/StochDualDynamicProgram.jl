# Risk Measures

A risk measure `ρ(X)` is a function that maps a set of random variables `X` to a real number.

The risk measures currently implemented are

- `Expectation()`

    ```julia
    SDDPModel(
        riskmeasure = Expectation()
    )
    ```

- `NestedCVaR(;lambda = 1, beta = 1)`

    `NestedCVaR[x] = lambda * E[x] + (1 - lambda) * CV@R(x, beta)`, where `CV@R(x, beta)` is the conditional value at risk of the beta quantile.

    For example, when minimising, `CV@R(x, 0.1)` is the expectation of the greatest 10% of outcomes.

    ```julia
    SDDPModel(
        riskmeasure = NestedCVaR(lambda=1, beta=0.5)
    )

    ```


## Creating a new risk measure

You can define a new risk measure by creating a new subtype of `AbstractRiskMeasure` and overloading a new `cutgenerator` method.

```julia
immutable MyRiskMeasure <: AbstractRiskMeasure
    # possibly containing fields
end

"""
    This function assembles a new cut using the following inputs
    + measure::AbstractRiskMeasure - used to dispatch
    + sense::Sense - either Maximum or Minimum
    + x::Vector{Vector{Float64}} - a vector of vector of state values for each scenario
    + pi::Vector{Vector{Float64}} - a vector of vector of dual values for each scenario
    + theta::Vector{Float64} - a vector of objective values for each scenario
    + prob::Vector{Float64} - the probability support of the scenarios. Should sum to one
    + stage::Int - the index of the stage
    + markov::Int - the index of the markov state
"""
function cutgenerator(measure::MyRiskMeasure, sense, x, pi, theta, prob, stage, markov)

end

```
