@everywhere begin
using Pkg; Pkg.activate("../")
end
using ArgParse
@everywhere using DrWatson
@everywhere begin
@quickactivate "EHTModCal" 
end
using Printf
@everywhere begin
    include("single_scan.jl")
end


@everywhere function parsim(dir, outdir, model; kwargs...)
    files = filter(endswith(".uvfits"), readdir(dir, join=true))
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
        nargs = '?'
        arg_type = Int
        default = 1
        help = "mring order to fit"
      "--dataset"
        nargs = '?'
        arg_type = Int
        default = 1
        help = "data set number"
      "--time"
        nargs = '?'
        arg_type = String
        default = "60"
        help = "Timescale to analyze"
    end

    parsed_args = parse_args(args,s)
    println(parsed_args)
    n = parsed_args["norder"]
    dataset = parsed_args["dataset"]
    stime = parsed_args["time"]

    if dataset < 35
      stage="uvfits_MoD_precal_stage1"
    elseif dataset < 119
      stage="uvfits_MoD_precal_stage2"
    else
      throw("Only stage 1 and 2 in right now $dataset !< 119")
    end
    dset = @sprintf "dataset%03d" dataset
    ddir = datadir("sims",
                   "P6ModCal",
                   stage,
                   dset,
                   "preimcal_pipeline",
                   stime,
                   "noisefrac0.02"
                  )
    outpath = projectdir("_research",
                        "P6ModCal_LCACP",
                        stage,
                        dset,
                        stime,
                        "noisefrac0.02"
                        )
    
    println("Order= $n")
    models = [ROSESoss.mring, ROSESoss.mringwgfloor]
    names  = ["mring","mring_gfloor"]
 
    for i in eachindex(models,names)
        parsim(ddir,
               joinpath(outpath,  names[i]*"_order-$n"),
               models[i](N=n,); 
               print_progress=false, dlogz=1.5
              )
    end
end


main(ARGS)




