using Pkg; Pkg.activate("../")
using ClusterManagers

using Distributed
using DrWatson
@quickactivate "EHTModCal" 


include("single_scan.jl")
include("launchruns.jl")



function massrun(orders, times, dsets)
  f = []
  for n in orders
        for t in times
          for d in dsets
            pids = addprocs(SlurmManager(40); nodes="1", partition="ehtq")

            @everywhere pids begin
              eval(:(using Pkg)) 
            end
            @everywhere pids Pkg.activate("../")
            @everywhere pids begin eval(:(using DrWatson)) end
            @everywhere pids begin
              @quickactivate "EHTModCal"
            end
            @everywhere pids include("single_scan.jl") 
            @everywhere pids include("launchruns.jl")
            task = remotecall(launchruns_val, pids[1], pids, d, n, t)
            push!(f, task)
          end
        end
      end
  
  return f
end



