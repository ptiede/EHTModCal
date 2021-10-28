using Pkg; Pkg.activate(".")
using ClusterManagers, Distributed
addprocs(SlurmManager(360), nodes="5", tasks_per_node="40", partition="ehtq")


@everywhere begin
  using Pkg; Pkg.activate(".")
end
@everywhere include(joinpath(@__DIR__, "stackermain.jl"))



function generatelist()
  dirs = filter(isdir, readdir("/gpfs/ptiede/EHTModCal/_research/GRMHDCal", join=true))
  dlist = String[]
  for d in dirs
      tmp = filter(x->occursin("gfloor", x), readdir(joinpath(d, "snapshot_120/noisefrac0.02"), join=true))
      push!(dlist, tmp...)
  end

  return dlist
end



function main()
  dlist = generatelist()
  pmap(dlist, on_error=e->println("Error $e")) do d
      if isfile(joinpath(d, "ChainHA_2", "chain_ha_trunc.csv"))
        println("Stacked chain exists skipping $d")
        return nothing
      end
      quickprocess(d)
  end
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end



