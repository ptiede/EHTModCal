using Pkg;Pkg.activate(joinpath(@__DIR__, "../../"))
include("plots.jl")
include("converters.jl")


function main()
    dirs = filter(startswith("../../GRMHDCal/hops_"), readdir("../../GRMHDCal/", join=true))
    paths = joinpath.(dirs, "snapshot_120/noisefrac0.02/mring_gfloor_order-4/ChainHA_HN/chain_ha_trunc.csv")
    outdir = mkpath("../../GRMHDCal/Summary/order4")
    summarize_ha.(paths, outdir)
end

main()
