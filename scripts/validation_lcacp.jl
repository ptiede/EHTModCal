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
    end

    parsed_args = parse_args(args,s)
    println(parsed_args)
    n = parsed_args["norder"]

    println("Order= $n")
    models = [ROSESoss.mringwfloor, ROSESoss.smringwfloor]
    names  = ["mring_floor", "mring_stretch_floor"]
    data = datadir("sims",
                   "ValidationData", 
                   "ehtim_thermal_preprocessing/")
    outpath =  projectdir("_research",
                      "ValidationRuns",
                      "ehtim_geometric_model",
                      "ehtim_thermal_preprocessing_lcacp") 
 
    for i in eachindex(models,names)
        parsim(data,
               joinpath(outpath,  names[i]*"_order-$n"),
               models[i](N=n,); print_progress=false, dlogz=1.5
              )
    end
end


main(ARGS)


