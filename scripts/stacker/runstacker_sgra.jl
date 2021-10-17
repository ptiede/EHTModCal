using Pkg; Pkg.activate(@__DIR__)
using Distributed, ClusterManagers
addprocs(SlurmManager(64), partition="amdq", nodes="2", tasks_per_node="32", t="0-23:59:59")

@everywhere begin
  using Pkg; Pkg.activate(".")
end
@everywhere include(joinpath(@__DIR__, "stackermain.jl"))


function generatelist()
    dlist = String[]
    for d in readdir("/gpfs/ptiede/EHTModCal/_research/SgrARunsFinal", join=true)
      for sd in readdir(d, join=true)
        for ssd in readdir(sd, join=true)
          for sssd in readdir(joinpath(ssd, "snapshot_120", "noisefrac0.02"), join=true)
            push!(dlist, sssd)
          end
        end
      end
    end
    dlist
end


function main()
  dlist = generatelist()
  pmap(dlist, on_error=e->println("Error $e")) do d
      if isfile(joinpath(d, "ChainHA_2", "chain_ha.csv"))
        println("Stacked chain exists skipping $d")
        return nothing
      end
      quickprocess(d)
  end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
