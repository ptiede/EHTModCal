#Activate correct environment
using Pkg; Pkg.activate(joinpath(@__DIR__,"."))
using Distributed
addprocs(4, exeflags="--project")
#This prevent Plots from trying to open a xserver window

#Do you want to plot the results?
@everywhere const plotr = false

@everywhere begin
if plotr
    ENV["GKSwstype"]="nul"
    using Plots
    using StatsPlots
end
end
#You only use dynesty right?
#using UltraNest


#Load ROSESoss
@everywhere begin
    using ROSESoss
    using DataFrames
    using CSV
    using DelimitedFiles
    using StatsBase
    using HDF5
    include(joinpath(@__DIR__,"fitscan.jl"))
end

"""
    makedirs(dir)
Makes the directory and subdirs that the scripts need
"""

@everywhere function makedirs(dir)
    mkpath(joinpath(dir))
    mkpath(joinpath(dir, "Chain"))
    mkpath(joinpath(dir, "MAP"))
    mkpath(joinpath(dir, "Residual"))
    mkpath(joinpath(dir, "Stats"))
end

@everywhere function singlescan(obsfile,
    outdir, scan;
    model = SossModels.smring2wfVACP,
    fitgains=true,
    kwargs...
    )
    makedirs(outdir)
    #Create output this does it relative to the directory this file is in
    open(joinpath(outdir, "params-$scan.dat"), "w") do io
    println(io, "Fitting Datafile: ", obsfile)
    println(io, "Using model: \n", model)
    println(io, "fit_gains = ", fitgains)
    end

    #Load data
    obs = ehtim.obsdata.load_uvfits(obsfile)
    obs.add_cphase(count="min")
    obs.add_amp(debias=true)
    #convert to ROSE EHTObservation objects
    ampobs = ROSESoss.extract_amps(obs)
    cpobs = ROSESoss.extract_cphase(obs)

    println("Start fit")
    dsum = fit_scan(outdir,
           model,
           ampobs, cpobs,
           scan,
           dynesty_sampler;
           fitgains = fitgains,
           plot_results=plotr,
           kwargs...,
           )
    println("Done fit")
    writedlm(joinpath(outdir, "Stats", "summary_stats-$scan.csv"), [keys(dsum), Tuple(dsum)], ',')
end

macro getmodel(m)
    sm = last(split("$m", "."))
    return sm
end


function parsim(dir, outdir, model; starti=1)
    files = filter(endswith(".uvfits"), readdir(dir, join=true))
    pmap(files) do f
        scan = parse(Int,first(split(last(split(f, "_snapshot")),".")))
        println("On scan $scan")
        singlescan(f, outdir, scan; model=model)
        nothing
    end
end


function mainrun()
    path = "synthetic_he/fullnight/ehtim_geometric_model/ehtim_thermal_phase_amp_preprocessing/"
    d = joinpath("ValidationData", path)
    parsim(d,
            joinpath("ValidationRuns", path,
            @getmodel(SossModels.smring1wfVACP)),
            SossModels.smring1wfVACP
          )
    parsim(d,
          joinpath("ValidationRuns", path,
          @getmodel(SossModels.smring2wfVACP)),
          SossModels.smring2wfVACP
        )
    parsim(d,
        joinpath("ValidationRuns", path,
        @getmodel(SossModels.smring3wfVACP)),
        SossModels.smring3wfVACP
      )
    parsim(d,
      joinpath("ValidationRuns", path,
      @getmodel(SossModels.smring4wfVACP)),
      SossModels.smring4wfVACP
    )
end
