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


function main()

    norder = 1:4
    models = [ROSESoss.mring, ROSESoss.mringwfloor, ROSESoss.smring, ROSESoss.smringwfloor]
    names  = ["mring", "mring_floor", "mring_stretch", "mring_stretch_floor"]
    data = datadir("sims",
                   "ValidationData", 
                   "ehtim_thermal_phase_amp_preprocessing",
                   "ehtim_thermal_phase_amp_preprocessing/")
    outpath =  projectdir("_research",
                      "ValidationRuns",
                      "ehtim_geometric_model",
                      "ehtim_thermal_phase_amp_preprocessing") 
 
    for i in eachindex(models,names)
      for n in norder
        parsim(data,
               joinpath(outpath,  names[i]*"_order-$n"),
               models[i](N=n,); print_progress=false
              )
      end
    end
end


main()

