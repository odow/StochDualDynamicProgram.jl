# Copyright Oscar Dowson, 2016

Cut(intercept, coefficients::Vector) = Cut(intercept, tuple(coefficients...))

calculatecut(riskmeasure, objective, values, duals) = error("Your risk measure needs to be an AbstractRiskMeasure. It's currently $(typeof(riskmeasure)).")

function calculatecut(riskmeasure::Expectation, objective, values, duals)
    Cut(objective - dot(duals, values), duals)
end

function _dot{N}(x::Tuple{Vararg{N, Float64}}, y::Tuple{Vararg{N, Float64}})
    z = 0.0
    @inbounds for i=1:N
        z += x[i] * y[i]
    end
    z
end
evaluate(cut::Cut{N}, state::Tuple{Vararg{N, Float64}}) = cut.intercept + _dot(cut.coefficients, state)

function addcut!{N}(sense::Sense, cs::CutStorage{N}, csinner::CutStorageInnter{N}, cut::Cut{N})
    push!(csinner.cuts, cut)
    push!(csinner.states_dominant, 0)
    for i in 1:length(cs.statesvisited)
        val = evaluate(cut, cs.statesvisited[i])
        if dominates(val, csinner.best_bound[i])
            csinner.best_bound[i] = val
            csinner.best_cut_index[i] = length(csinner.cuts)
            csinner.states_dominant[end] += 1
        end
    end
end

# function rebuild!(::DeMatosCutSelection, m::SDDPModel)
# end

# TODO: need to be able to run cut selection on a single rib
