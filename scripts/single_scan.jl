#Activate correct environment
#This prevent Plots from trying to open a xserver window


#Load ROSESoss
using DrWatson
@quickactivate "."
using ROSESoss
using DataFrames
using CSV
using DelimitedFiles
using StatsBase
using HDF5
using HypercubeTransform
using NamedTupleTools
using TupleVectors
include(scriptsdir("fitscan.jl"))


"""
    makedirs(dir)
Makes the directory and subdirs that the scripts need
"""
function makedirs(dir)
    mkpath(joinpath(dir))
    mkpath(joinpath(dir, "Chain"))
    mkpath(joinpath(dir, "MAP"))
    mkpath(joinpath(dir, "Residual"))
    mkpath(joinpath(dir, "Stats"))
end



function singlescan_vacp(obsfile,
                 outdir, scan,
                 model;
                 tmin=0.0, 
                 tmax=24.0,
                 kwargs...
                 )

    makedirs(outdir)

    isfile(joinpath(outdir, "Chain",  "nested_chain_scan-$scan.csv")) && return nothing

    #Load data
    obs = ehtim.obsdata.load_uvfits(obsfile)
    obs.add_amp(debias=true)
    #convert to ROSE EHTObservation objects
    ampobs = ROSESoss.extract_amps(obs)
    if (ampobs[1].time < tmin) || (ampobs[1].time > tmax)
      println("Not in good time skipping")
      return nothing
    end
    bl = ROSE.getdata(ampobs,:baselines)
    s1 = unique(first.(bl))
    s2 = unique(last.(bl))
    stations = unique([s1...,s2...])
    if length(stations) < 4
      return nothing
    end

    obs.add_cphase(count="min-cut0bl")
    cpobs = ROSESoss.extract_cphase(obs)

    #Create output this does it relative to the directory this file is in
    open(joinpath(outdir, "params-$scan.dat"), "w") do io
        println(io, "Fitting Datafile: ", obsfile)
        println(io, "Using model: \n", model)
        println(io, "stations: ", stations...)
    end


    println("Start fit")
    chain, echain, stats =
                fit_scan(outdir,
                    model,
                    ampobs, cpobs,
                    DynestyStatic(nlive=2000);
                    kwargs...,
                    )
    println("Done fit writing results to disk")
    CSV.write(joinpath(outdir, "Stats", "summary_stats-$scan.csv"), nt2df(stats))
    #tv2df(chain) |> CSV.write(joinpath(outdir, "Chain",  "nested_chain_scan-$scan.csv"))
    tv2df(echain) |> CSV.write(joinpath(outdir, "Chain",  "equal_chain_scan-$scan.csv"))

    nothing
end


function singlescan_lcacp(obsfile,
    outdir, scan,
    model;
    kwargs...
    )
    makedirs(outdir)
    println(outdir)
    #Create output this does it relative to the directory this file is in
    open(joinpath(outdir, "params-$scan.dat"), "w") do io
      println(io, "Fitting Datafile: ", obsfile)
      println(io, "Using model: \n", model)
    end

    #Load data
    obs = ehtim.obsdata.load_uvfits(obsfile)

    try
      obs.add_logcamp(debias=true, count="min")
    catch
      @warn "Warning no quadrangle so skipping scan $scan"
      return nothing
    end
     obs.add_cphase(count="min-cut0bl")




    ampobs = ROSESoss.extract_lcamp(obs)
    cpobs = ROSESoss.extract_cphase(obs)


    open(joinpath(outdir, "params-$scan.dat"), "w") do io
        println(io, "Fitting Datafile: ", obsfile)
        println(io, "Using model: \n", model)
    end

    println("Start fit")
    chain, echain, stats =
                fit_scan(outdir,
                    model,
                    ampobs, cpobs,
                    DynestyStatic(nlive=2000);
                    kwargs...,
                    )
    println("Done fit writing results to disk")
    CSV.write(joinpath(outdir, "Stats", "summary_stats-$scan.csv"), nt2df(stats))
    #tv2df(chain) |> CSV.write(joinpath(outdir, "Chain",  "nested_chain_scan-$scan.csv"))
    tv2df(echain) |> CSV.write(joinpath(outdir, "Chain",  "equal_chain_scan-$scan.csv"))
    nothing
end


macro getmodel(m)
    sm = quote last(split("$(m)", ".")) end
    return sm
end
