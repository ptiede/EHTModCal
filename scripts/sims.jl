@everywhere begin
using Pkg; Pkg.activate("../")
end
@everywhere using DrWatson
@everywhere begin
@quickactivate "EHTModCal" 
end

@everywhere begin
    include("single_scan.jl")
end

@everywhere function sim(dir, outdir, model; kwargs...)
    files = filter(endswith(".uvfits"), readdir(dir, join=true))
    for f in files
        scan = parse(Int,first(split(last(split(basename(f), "_snapshot")),".")))
        println("On scan $scan")
        singlescan_vacp(f, outdir, scan, model; kwargs...)
    end
end

@everywhere function parsim(dir, outdir, model; kwargs...)
    files = filter(endswith(".uvfits"), readdir(dir, join=true))
    pmap(files) do f
        scan = parse(Int,first(split(last(split(basename(f), "_snapshot")),".")))
        println("On scan $scan")
        singlescan_vacp(f, outdir, scan, model; kwargs...)
        GC.gc()
    end
end
