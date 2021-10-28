using Pkg
Pkg.activate("/scratch/ptiede/SgrA/EHTModCal/scripts/stacker")

@everywhere begin
  using Pkg; Pkg.activate("/scratch/ptiede/SgrA/EHTModCal/scripts/stacker")
end
@everywhere include(joinpath("/scratch/ptiede/SgrA/EHTModCal/scripts/stacker", "stackermain.jl"))

using Comonicon



function generatelist(order)
  dirs = filter(x->startswith(basename(x), "hops"), readdir("/scratch/ptiede/SgrA/EHTModCal/GRMHDVal", join=true))
  dlist = String[]
  for d in dirs
      tmp = filter(x->occursin("gfloor_order-$(order)", x), readdir(joinpath(d, "snapshot_120/noisefrac0.02"), join=true))
      push!(dlist, tmp...)
  end

  return dlist
end


"""
Runs the GRMHD stacker for some m-ring order

# Arguments

- `o`: The mring order you can to analyze
"""
@main function main(o::Int)
  dlist = generatelist(o)
  println(dlist)
  pmap(dlist) do d
      quickprocess(d)
  end
end
