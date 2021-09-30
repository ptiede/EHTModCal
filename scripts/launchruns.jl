using Printf

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

function launchruns_sgra(pids, order, data, time, day, band; kwargs...)
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
       parsim(pid[2:end], ddir,
               joinpath(outpath,  names[i]*"_order-$n"),
               models[i](N=n,); 
               print_progress=false, dlogz=1.5
              )
    end
  
end


function launchruns_val(pids, dataset::Int, order, time; kwargs...)
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
                   time,
                   "noisefrac0.02"
                  )
    outpath = projectdir("_research",
                        "P6ModCal",
                        stage,
                        dset,
                        time,
                        "noisefrac0.02"
                        )
    
    n = order
    println("Order= $n")
    models = [ROSESoss.mring, ROSESoss.mringwfloor, ROSESoss.smringwfloor]
    names  = ["mring", "mring_floor", "mring_floor_stretch"]
    
    println("starting runs and outputting to $(outpath)")
    for i in eachindex(models,names)
      parsim(pids[2:end], ddir,
               joinpath(outpath,  names[i]*"_order-$n"),
               models[i](N=n,); 
               print_progress=false, dlogz=1.5
              )
    end
    #for p in pids[2:end]
    #    try
    #    remotecall_wait(exit, p)
    #  catch
    #    println("exiting $p")
    #  end
    #end
    #try
    #  remotecall(exit, pids[1])
    #catch
    #  println("We are done")
    #end
end
