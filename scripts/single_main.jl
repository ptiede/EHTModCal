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
    mkpath(joinpath(@__DIR__, dir))
    mkpath(joinpath(@__DIR__, dir, "Chain"))
    mkpath(joinpath(@__DIR__, dir, "MAP"))
    mkpath(joinpath(@__DIR__, dir, "Residual"))
    mkpath(joinpath(@__DIR__, dir, "Stats"))
end

function multiscanvis(obsfile, outdir; model=SossModels.smring1wfVis, fitgains=true)
    makedirs(outdir)
    #Create output this does it relative to the directory this file is in
    open(joinpath(@__DIR__, outdir, "params.dat"), "w") do io
        println(io, "Fitting Datafile: ", obsfile)
        println(io, "Using model: \n", model)
        println(io, "fit_gains = ", fitgains)
    end

    obs = ehtim.obsdata.load_uvfits(obsfile)
    obs_scans = obs.split_obs()
    println("There are $(length(obs_scans)) scans")
    vis = extract_vis.(obs_scans)

    map(1:length(vis)) do i
        dsum = fit_scan(outdir,
                 model,
                 vis[i],
                 i-1,
                 dynesty_sampler;
                 fitgains = fitgains,
                 plot_results=plotr,
                )
        println("Done fit")
        writedlm(joinpath(@__DIR__, outdir, "Stats", "summary_stats-$(i-1).csv"), [keys(dsum), Tuple(dsum)], ',')
    end
end



function multiscan(obsfile, outdir; model=SossModels.smring2wfVACP, fitgains=true)
    makedirs(outdir)
    #Create output this does it relative to the directory this file is in
    open(joinpath(@__DIR__, outdir, "params.dat"), "w") do io
        println(io, "Fitting Datafile: ", obsfile)
        println(io, "Using model: \n", model)
        println(io, "fit_gains = ", fitgains)
    end

    obs = ehtim.obsdata.load_uvfits(obsfile)
    obs_scans = obs.split_obs()
    for o in obs_scans
        o.add_amp(debias=true)
        o.add_cphase(count="min")
    end
    println("There are $(length(obs_scans)) scans")
    amps = extract_amps.(obs_scans)
    phases = extract_cphase.(obs_scans)

    pmap(1:length(amps)) do i
        dsum = fit_scan(outdir,
                 model,
                 amps[i],
                 phases[i],
                 i-1,
                 dynesty_sampler;
                 fitgains = fitgains,
                 plot_results=plotr,
                )
        println("Done fit")
        writedlm(joinpath(@__DIR__, outdir, "Stats", "summary_stats-$(i-1).csv"), [keys(dsum), Tuple(dsum)], ',')
    end



end

@everywhere function singlescan(obsfile,
                 outdir, scan;
                 model = SossModels.smring2wfVACP,
                 fitgains=true,
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
    writedlm(joinpath(@__DIR__, outdir, "Stats", "summary_stats-$scan.csv"), [keys(dsum), Tuple(dsum)], ',')
end





@everywhere function sim(dir, outdir, model)
    files = readdir(dir, join=true)
    for f in files[134:end]
        scan = parse(Int,first(split(last(split(f, "_scan")),".")))
        println("On scan $scan")
        singlescan(f, outdir, scan; model=model)
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

function modmain(m = SossModels.smring3wfVACP)
    dirs = filter(isdir, readdir("MoD_calibration_suite/", join=true))
    @sync @distributed for d in dirs
        println("On directory $d")
        sub = "preprocessed_test/scan-average/noisefrac0.02/"
        od = joinpath(@__DIR__, "Runs", last(splitpath(d)))
        sim(joinpath(d, sub), joinpath(od,sub), m)
    end
end

function fullmain()
    d = "SgrAData/SgrA_ER6_ver2021-04-15/3599/hops/hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s/preimcal_pipeline_camp/"
    o = "RunsCal/ROSERuns/SgrARuns/3599"
    ts = ["scan-average"]
    for t in ts
        parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.mring1wfVACP)), SossModels.mring1wfVACP)
        parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.smring1wfVACP)), SossModels.smring1wfVACP)
        parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.mring2wfVACP)), SossModels.mring2wfVACP)
        parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.smring2wfVACP)), SossModels.smring2wfVACP)
        #parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.mring3wfVACP)), SossModels.mring3wfVACP)
        #parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.smring3wfVACP)), SossModels.smring3wfVACP)
    end

    #d = "SgrAData/SgrA_ER6_ver2021-04-15/3598/hops/hops_3598_SGRA_LO_netcal_LMTcal_normalized_10s/preimcal_pipeline_camp/"
    #o = "RunsCal/ROSERuns/SgrARuns/3598"
    #ts = ["scan-average", "snapshot60","snapshot120"]
    #for t in ts
        #parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.mring2wfVACP)), SossModels.mring2wfVACP)
        #parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.smring2wfVACP)), SossModels.smring2wfVACP)
        #parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.mring3wfVACP)), SossModels.mring3wfVACP)
        #parsim(joinpath(d,t,"noisefrac0.02"), joinpath(o, t, @getmodel(SossModels.smring3wfVACP)), SossModels.smring3wfVACP)
    #end
end


function bt98main()
    d98 = "SgrAData/BestTimes/hops_3598_SGRA_LO_netcal_LMTcal_10s_regionIII/preimcal_pipeline_camp/snapshot60/noisefrac0.02/"
    o98 = "RunsCal/SgrARuns/BestTimes/3598/preimcal_camp/snapshot60/"

    parsim(d98, joinpath(o98, @getmodel(SossModels.mring1wfVACP)), SossModels.mring1wfVACP)
    parsim(d98, joinpath(o98, @getmodel(SossModels.mring2wfVACP)), SossModels.mring2wfVACP)
    parsim(d98, joinpath(o98, @getmodel(SossModels.mring3wfVACP)), SossModels.mring3wfVACP)
    parsim(d98, joinpath(o98, @getmodel(SossModels.mring4wfVACP)), SossModels.mring4wfVACP)
    parsim(d98, joinpath(o98, @getmodel(SossModels.smring1wfVACP)), SossModels.smring1wfVACP)
    parsim(d98, joinpath(o98, @getmodel(SossModels.smring2wfVACP)), SossModels.smring2wfVACP)
    parsim(d98, joinpath(o98, @getmodel(SossModels.smring3wfVACP)), SossModels.smring3wfVACP)
    parsim(d98, joinpath(o98, @getmodel(SossModels.smring4wfVACP)), SossModels.smring4wfVACP)
end

function bt99main()
    d99 = "SgrAData/BestTimes/hops_3599_SGRA_LO_netcal_LMTcal_10s_regionIII/preimcal_pipeline_camp/snapshot60/noisefrac0.02/"
    o99 = "RunsCal/SgrARuns/BestTimes/3599min/preimcal_camp/snapshot60/"

    #parsim(d99, joinpath(o99, @getmodel(SossModels.mring1wfVACP)), SossModels.mring1wfVACP)
    #parsim(d99, joinpath(o99, @getmodel(SossModels.mring2wfVACP)), SossModels.mring2wfVACP)
    parsim(d99, joinpath(o99, @getmodel(SossModels.mring3wfVACP)), SossModels.mring3wfVACP)
    #parsim(d99, joinpath(o99, @getmodel(SossModels.mring4wfVACP)), SossModels.mring4wfVACP)
    #parsim(d99, joinpath(o99, @getmodel(SossModels.smring1wfVACP)), SossModels.smring1wfVACP)
    #parsim(d99, joinpath(o99, @getmodel(SossModels.smring2wfVACP)), SossModels.smring2wfVACP)
    #parsim(d99, joinpath(o99, @getmodel(SossModels.smring3wfVACP)), SossModels.smring3wfVACP)
    #parsim(d99, joinpath(o99, @getmodel(SossModels.smring4wfVACP)), SossModels.smring4wfVACP)
end


function overnight()
    dirfol = "SgrAData/SgrA_ER6_ver2021-04-15/3599/hops/hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s/preimcal_pipeline_camp/"


    parsim(joinpath(dirfol, "snapshot120/noisefrac0.02/"),
           "RunsCal/SgrARuns/preimcal_camp/snapshot120/hops/Mring3wfVACP",
           SossModels.mring3wfVACP, starti=134)

    parsim(joinpath(dirfol, "snapshot60/noisefrac0.02/"),
           "RunsCal/SgrARuns/preimcal_camp/snapshot60/hops/Mring3wfVACP",
           SossModels.mring3wfVACP)



end
