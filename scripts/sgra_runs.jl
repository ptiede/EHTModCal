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


@everywhere function parsim(dir, outdir, model; kwargs...)
    files = filter(endswith(".uvfits"), readdir(dir, join=true))
    pmap(files, on_error=e->(println("Error on $dir"))) do f
        scan = parse(Int,first(split(last(split(basename(f), "_scan")),".")))
        println("On scan $scan")
        singlescan_vacp(f, outdir, scan, model; kwargs...)
        GC.gc()
    end
end


function main(args)

    s = ArgParseSettings(description="Run the Validation set with different mring orders")

    @add_arg_table! s begin
      "--norder"
        nargs = '?'
        arg_type = Int
        default = 1
        help = "mring order to fit"
      "--data"
        nargs = '?'
        arg_type = String
        default = "hops"
        help = "data set either hops or casa"
      "--time"
        nargs = '?'
        arg_type = String
        default = "60"
        help = "Timescale to analyze"
      "--day"
        nargs = '?'
        arg_type = String
        default = "3599"
        help = "observation day to analyze"
       "--band"
        nargs = '?'
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
    outpath =  DrWatson.projectdir("_research",
                      "SgrARunsFinal",
                      data,
                      day,
                      band,
                      "snapshot_$stime",
                      "noisefrac0.02")
    println("Saving to $(outpath)")
 
    for i in eachindex(models,names)
        parsim(ddir,
               joinpath(outpath,  names[i]*"_order-$n"),
               models[i](N=n,); 
               print_progress=false
              )
    end
end


main(ARGS)


