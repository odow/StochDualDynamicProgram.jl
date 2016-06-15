const filedir  = dirname(@__FILE__)
const jsonfile = joinpath(filedir, "run.json")
const htmlfile = joinpath(filedir, "Visualise.html")

function addtoputput!(output::Dict{ASCIIString, Any}, sym, value)
	for key in keys(output)
		if string(sym) == key
			output[key] = value
			return
		end
	end
	error("Keyword $(sym)=$value not recognised in @visualise.")
end

macro visualise(results, kw, block)
	@assert block.head == :block || error("Invalid syntax for @visualise")
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
			output = Dict{ASCIIString, Any}(
			    "cumulative"  => false,
			    "title"       => "",
			    "ylabel"      => "",
			    "xlabel"      => "Stages",
			    "interpolate"  => "linear"
			)
            if it.head == :tuple
                if length(it.args) > 2
                    error("Unknown arguments in @visualise")
                end

				if length(it.args) == 2
					if it.args[2].head == :tuple
						for arg in it.args[2].args
							if arg.head != :(=)
								error("Must be a keyword argument in @visualise: $(arg)")
							end
							addtoputput!(output, arg.args[1], arg.args[2])
						end
					elseif it.args[2].head == :(=)
						addtoputput!(output, it.args[2].args[1],it.args[2].args[2])
					end
					f = Expr(:->, kw, Expr(:block, it.args[1]))
				else
					f = Expr(:->, kw, Expr(:block, it.args))
				end
            else
				f = Expr(:->, kw, Expr(:block, it))
            end
			push!(code.args, quote
				adddata!(visualiseout, $(esc(results)), replications, stages, $(esc(f)), $output)
			end)
        end
    end
	tmphtmlfile = replace(tempname(), ".tmp", ".html")
	runcmd = `$(ENV["COMSPEC"]) /c start $tmphtmlfile`
	visualisehtml = readall(htmlfile)
	push!(code.args, quote
		if OS_NAME == :Windows
			open($tmphtmlfile, "w") do f
				write(f, replace($visualisehtml, "<!--DATA-->", json(visualiseout)))
			end
			run($runcmd)
		else
			replace($visualisehtml, "<!--DATA-->", json(visualiseout))
		end
	end)
    return code
end

function adddata!(visualiseout::Vector{Dict{ASCIIString, Any}}, results::Dict{Symbol, Any}, replications::Int, stages::Int, func::Function, output::Dict{ASCIIString, Any})
	output["data"] = Array(Vector{Float64}, replications)
	for i=1:replications
		output["data"][i] = zeros(stages)
		y = 0.
		for j=1:stages
			if output["cumulative"]
				y += func(j, i, results)
			else
				y = func(j, i, results)
			end
			output["data"][i][j] = y
		end
	end
	push!(visualiseout, output)
end
