using Pkg; Pkg.activate(joinpath(@__DIR__, "../"))
using ClusterManagers, Distributed
addprocs(SlurmManager(200), nodes="5", tasks_per_node="40", partition="ehtq")
include(joinpath(@__DIR__,"create_grmhddirs.jl"))

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


@everywhere function parsim(files, outdir, model; kwargs...)
  pmap(zip(files, outdir)) do (f,o)
        scan = parse(Int,first(split(last(split(basename(f), "_scan")),".")))
        println("On scan $scan")
        singlescan_vacp(f, o, scan, model; kwargs...)
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
    day = parsed_args["day"]
    band = parsed_args["band"]

    println("Order= $n")
    models = [ROSESoss.mring, ROSESoss.mringwgfloor]
    names  = ["mring", "mring_gfloor"]

    files = filter(endswith(".uvfits"), create_filelist(create_paths(band, day)))
    odirs = create_outdirs(files)
    println("starting fits")

    for i in eachindex(models,names)
        println("Fitting model $(names[i])")
        parsim(files,
               joinpath.(odirs,  Ref(names[i]*"_order-$n")),
               models[i](N=n, fmin=0.8, fmax=1.2);
               print_progress=false
              )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end

