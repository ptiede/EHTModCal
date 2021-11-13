using Pkg
Pkg.activate("/scratch/ptiede/SgrA/EHTModCal/scripts/stacker")

@everywhere begin
  using Pkg; Pkg.activate("/scratch/ptiede/SgrA/EHTModCal/scripts/stacker")
end
@everywhere include(joinpath("/scratch/ptiede/SgrA/EHTModCal/scripts/stacker", "stackermain.jl"))

using Comonicon





function generatelist(order)
    dirs = [
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_MAD_a-0.5_Rh160_i10_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_MAD_a-0.94_Rh10_i50_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_MAD_a+0.5_Rh40_i50_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_MAD_a+0.5_Rh160_i90_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_MAD_a+0.94_Rh10_i90_LO/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_MAD_a0_Rh160_i10_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a-0.5_Rh40_i10_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a-0.5_Rh160_i90_LO/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a-0.94_Rh40_i90_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a+0.5_Rh10_i10_LO/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a+0.94_Rh160_i90_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a0_Rh10_i10_LO/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a0_Rh10_i50_HI/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a0_Rh10_i90_LO/",
        "/scratch/ptiede/SgrA/EHTModCal/GRMHDCal/hops_3599_SANE_a0_Rh40_i50_LO/",
    ]
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
      quickprocess(d; restart=true, nsteps=5_000_000)
  end
end
