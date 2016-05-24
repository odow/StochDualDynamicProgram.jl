using JSON

const filedir  = dirname(@__FILE__)
const jsonfile = joinpath(filedir, "run.json")
const htmlfile = joinpath(filedir, "Visualise.html")

macro visualise(results, kw, block)
	@assert block.head == :block || error("Invalid syntax for @visualisesimulation")
	push!(kw.args, results)
    code = quote
		visualiseout = Dict{ASCIIString, Any}[]
		replications = length($(esc(results))[:Current][1])
		stages = length($(esc(results))[:Current])
	end
    for it in block.args
        if Base.Meta.isexpr(it, :line)
            # do nothing
        else
            if it.head == :tuple
				iscumulative = false
                if length(it.args) == 3
					if it.args[3] == :(cumulative = true)
						iscumulative = true
					end
				elseif length(it.args) > 3
                    error("Unknown arguments in @visualisesimulation")
                end
				f = Expr(:->, kw, Expr(:block, it.args[1]))
                push!(code.args, quote
                    adddata!(visualiseout, $(esc(results)), replications, stages, $(esc(f)), $iscumulative, $(it.args[2]))
                end)
            else
                error("Unknown arguments in @visualisesimulation")
            end
        end
    end
	tmphtmlfile = replace(tempname(), ".tmp", ".html")
	runcmd = `$(ENV["COMSPEC"]) /c start $tmphtmlfile`
	visualisehtml = readall(htmlfile)
	push!(code.args, quote
		open($tmphtmlfile, "w") do f
			write(f, replace($visualisehtml, "<!--DATA-->", json(visualiseout)))
		end
		run($runcmd)
	end)
    return code
end

function adddata!(visualiseout::Vector{Dict{ASCIIString, Any}}, results::Dict{Symbol, Any}, replications::Int, stages::Int, func::Function, iscumulative::Bool, label)
	output = Dict{ASCIIString, Any}()
	output["label"] = label
	output["data"] = Array(Vector{Float64}, replications)
	for i=1:replications
		output["data"][i] = zeros(stages)
		y = 0.
		for j=1:stages
			if iscumulative
				y += func(j, i, results)
			else
				y = func(j, i, results)
			end
			output["data"][i][j] = y
		end
	end
	push!(visualiseout, output)
end
