#Activate correct environment
using Pkg; Pkg.activate(joinpath(@__DIR__,"."))
using Distributed
addprocs(2, exeflags="--project")
#This prevent Plots from trying to open a xserver window


#Load ROSESoss
@everywhere begin
    using ROSESoss
    using DataFrames
    using CSV
    using DelimitedFiles
    using StatsBase
    using HDF5
    using HyercubeTransform
    using NamedTupleTools
    include(joinpath(@__DIR__,"fitscan.jl"))

end

@everywhere begin
"""
    makedirs(dir)
Makes the directory and subdirs that the scripts need
"""
function makedirs(dir)
    mkpath(joinpath(@__DIR__, dir))
    mkpath(joinpath(@__DIR__, dir, "Chain"))
    mkpath(joinpath(@__DIR__, dir, "MAP"))
    mkpath(joinpath(@__DIR__, dir, "Residual"))
    mkpath(joinpath(@__DIR__, dir, "Stats"))
end

end #@eveywhere


@everywhere function singlescan_vacp(obsfile,
                 outdir, scan,
                 model;
                 kwargs...
                 )
    makedirs(outdir)
    #Create output this does it relative to the directory this file is in
    open(joinpath(@__DIR__, outdir, "params-$scan.dat"), "w") do io
        println(io, "Fitting Datafile: ", obsfile)
        println(io, "Using model: \n", model)
        println(io, "fit_gains = ", fitgains)
    end

    #Load data
    obs = ehtim.obsdata.load_uvfits(obsfile)
    obs.add_cphase(count="min-cut0bl")
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
    writedlm(joinpath(@__DIR__, outdir, "Stats", "summary_stats-$scan.csv"), [keys(dsum), Tuple(dsum)], ',')
end


@everywhere function singlescan_lcacp(obsfile,
    outdir, scan,
    model;
    kwargs...
    )
    makedirs(outdir)
    #Create output this does it relative to the directory this file is in
    open(joinpath(@__DIR__, outdir, "params-$scan.dat"), "w") do io
    println(io, "Fitting Datafile: ", obsfile)
    println(io, "Using model: \n", model)
    println(io, "fit_gains = ", fitgains)
    end

    #Load data
    obs = ehtim.obsdata.load_uvfits(obsfile)
    obs.add_cphase(count="min-cut0bl")
    obs.add_logcamp(debias=true, count="min")
    #convert to ROSE EHTObservation objects
    ampobs = ROSESoss.extract_lcamps(obs)
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
    writedlm(joinpath(@__DIR__, outdir, "Stats", "summary_stats-$scan.csv"), [keys(dsum), Tuple(dsum)], ',')
end





@everywhere function sim(dir, outdir, model)
    files = readdir(dir, join=true)
    for f in files[134:end]
        scan = parse(Int,first(split(last(split(f, "_scan")),".")))
        println("On scan $scan")
        singlescan_vacp(f, outdir, scan, model)
    end
end

macro getmodel(m)
    sm = last(split("$m", "."))
    return sm
end


function parsim(dir, outdir, model; starti=1)
    files = readdir(dir, join=true)
    pmap(files) do f
        scan = parse(Int,first(split(last(split(f, "_scan")),".")))
        println("On scan $scan")
        singlescan(f, outdir, scan; model=model)
        nothing
    end
end
