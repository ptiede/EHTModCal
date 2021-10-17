@everywhere begin
  using Pkg; Pkg.activate(joinpath(@__DIR__, "../"))
end
using ArgParse
@everywhere using DrWatson
@everywhere begin
@quickactivate "EHTModCal"
end

@everywhere begin
    include("single_scan.jl")
end


function parsim(pids, dir, outdir, model; kwargs...)
    files = filter(endswith(".uvfits"), readdir(dir, join=true))
    pmap(WorkerPool(pids), files) do f
        scan = parse(Int,first(split(last(split(basename(f), "_scan")),".")))
        println("On scan $scan")
        singlescan_vacp(f, outdir, scan, model; kwargs...)
        GC.gc()
    end
end

function parsimorder(pids, dir, outdir, orders; kwargs...)
    for n in orders
      outfile = joinpath(outdir, "mring_gfloor_order-$n")
      m = ROSESoss.mringwgfloor(N=3, fmin=0.1, fmax=4.0)
      parsim(pids, dir, outfile, m; kwargs...)
    end
end



function parsimfiles(pids, files, outdir, model; kwargs...)
  pmap(WorkerPool(pids), files) do f
        scan = parse(Int,first(split(last(split(basename(f), "_scan")),".")))
        println("On scan $scan")
        singlescan_lcacp(f, outdir, scan, model; kwargs...)
        GC.gc()
    end
end
