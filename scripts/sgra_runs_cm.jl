using Pkg; Pkg.activate("../")
using ClusterManagers

using Distributed
using DrWatson
@quickactivate "EHTModCal" 


include("single_scan.jl")


function parsim(pids, dir, outdir, model; kwargs...)
    files = filter(endswith(".uvfits"), readdir(dir, join=true))
    println("start pmap")
    pool = WorkerPool(pids)
    println(pool)
    pmap(WorkerPool(pids), files) do f
        scan = parse(Int,first(split(last(split(basename(f), "_scan")),".")))
        println("On scan $scan")
        singlescan_vacp(f, outdir, scan, model; kwargs...)
        GC.gc()
    end
end


function launchruns(nworkers, order, data, time, day, band; kwargs...)
    pids = addprocs(SlurmManager(nworkers); nodes="1", kwargs...)
    @everywhere pids begin
      eval(:(using Pkg)) 
    end
    @everywhere pids Pkg.activate("../")
    @everywhere pids eval(:(using DrWatson))
    @everywhere pids begin
      @quickactivate "EHTModCal"
    end
    @everywhere pids include("single_scan.jl") 
    n = order
    println("Order= $n")
    models = [ROSESoss.mring, ROSESoss.mringwfloor, ROSESoss.smringwfloor]
    names  = ["mring", "mring_floor", "mring_floor_stretch"]
    ddir = datadir("exp_pro",
                   "preprocessed_data", 
                   data,
                   day,
                   band,
                   "snapshot_$time",
                   "noisefrac0.02")
    outpath =  projectdir("_research",
                      "SgrARuns",
                      data,
                      day,
                      band,
                      "snapshot_$time",
                      "noisefrac0.02") 
    println("starting runs and outputting to $(outpath)")
    for i in eachindex(models,names)
        parsim(pids, ddir,
               joinpath(outpath,  names[i]*"_order-$n"),
               models[i](N=n,); 
               print_progress=false, dlogz=1.5
              )
    end
    rmprocs(pids)
end


function letsgo(orders, times, bands, days, pipeline)

  f = @async begin
    for d in days
      @async for b in bands
        @async for t in times
          @async for n in orders
            @async launchruns(32, n, pipeline, t, d, b, partition="ehtq")
          end
        end
      end
    end
  end
  return f
end



