using Pkg;Pkg.activate(joinpath(@__DIR__, "../../"))
include("plots.jl")
include("converters.jl")


function main()
    dirs = filter(startswith("../../GRMHDVal/hops_"), readdir("../../GRMHDVal/", join=true))
    paths = joinpath.(dirs, "snapshot_120/noisefrac0.02/mring_gfloor_order-3/ChainHA_2/chain_ha_trunc.csv")
    outdir = mkpath("../../GRMHDVal/Summary/order3")
    summarize_ha.(paths, outdir)

    paths = joinpath.(dirs, "snapshot_120/noisefrac0.02/mring_gfloor_order-4/ChainHA_2/chain_ha_trunc.csv")
    outdir = mkpath("../../GRMHDVal/Summary/order4")
    summarize_ha.(paths, outdir)
end

#main()
