@everywhere begin
using Pkg; Pkg.activate("../")
end
using ArgParse
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
        scan = parse(Int,first(split(last(split(basename(f), "_scan")),".")))
        println("On scan $scan")
        singlescan_lcacp(f, outdir, scan, model; kwargs...)
    end
end

@everywhere function parsim(dir, outdir, model; kwargs...)
    files = filter(endswith(".uvfits"), readdir(dir, join=true))
    println("There are $(length(files))")
    pmap(files) do f
        scan = parse(Int,first(split(last(split(basename(f), "_scan")),".")))
        println("On scan $scan")
        singlescan_lcacp(f, outdir, scan, model; kwargs...)
        GC.gc()
    end
end


function main(args)

    s = ArgParseSettings(description="Run the Validation set with different mring orders")

    @add_arg_table! s begin
      "--norder"
        arg_type = Int
        default = 1
        help = "mring order to fit"
      "--data"
        arg_type = String
        default = "hops"
        help = "data set either hops or casa"
      "--time"
        arg_type = String
        default = "60"
        help = "Timescale to analyze"
      "--day"
        arg_type = String
        default = "3599"
        help = "observation day to analyze"
       "--band"
        arg_type = String
        default = "LO"
        help = "observation band to use"
    end

    parsed_args = parse_args(args,s)
    println(parsed_args)
    n = parsed_args["norder"]
    data = parsed_args["data"]
    stime = parsed_args["time"]
    day = parsed_args["day"]
    band = parsed_args["band"]

    println("Order= $n")
    
    models = [ROSESoss.mring, ROSESoss.mringwgfloor]
    names  = ["mring", "mring_gfloor"]
    ddir = datadir("exp_pro",
                   "preprocessed_data", 
                   data,
                   day,
                   band,
                   "snapshot_$stime",
                  "noisefrac0.02")
    outpath =  projectdir("_research",
                      "SgrARunsFinal_LCACP",
                      data,
                      day,
                      band,
                      "snapshot_$stime") 
 
    for i in eachindex(models,names)
        parsim(ddir,
               joinpath(outpath,  names[i]*"_order-$n"),
               models[i](N=n,); print_progress=false, dlogz=1.5
              )
    end
    
end


main(ARGS)


