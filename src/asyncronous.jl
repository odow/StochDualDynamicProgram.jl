function asyncsolve()
    np = nprocs()  # determine the number of processes available

    cuts = Cut[]
    indices = zeros(np)

    converged = false
    sddpisconverged() = (isconverged=converged; isconverged)
    setconverged!() = (converged=true)
    @sync begin
        for p=1:np # for each process
            if p != myid() || np == 1 # exclude the master process
                @async begin    # async
                    while !sddpisconverged()
                        (discoveredcuts, new_convergence) = remotecall_fetch(getnewcuts!, p, cuts[indices[p]:end])
                        for cut in discoveredcuts
                            push!(cuts, cut)
                        end
                        indices[p] = length(cuts)
                        if new_convergence
                            setconverged!()
                        end
                    end
                end
            end
        end
    end
end
