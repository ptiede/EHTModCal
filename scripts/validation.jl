using Pkg; Pkg.activate("")
using DrWatson
@quickactivate "."
using Distributed
addprocs(2, exeflags="--project")

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
        singlescan_vacp(f, outdir, scan, model; kwargs...)
    end
end


function main()
    parsim(datadir("sims","ValidationData/ehtim_thermal_phase_amp_preprocessing/ehtim_thermal_phase_amp_preprocessing/"),
           projectdir("_research","ValidationRuns","ehtim_geometric_model","ehtim_thermal_phase_amp_preprocessing"),
           ROSESoss.mring(N=3,))
end



function mainrun()
    path = "synthetic_he/fullnight/ehtim_geometric_model/ehtim_thermal_phase_amp_preprocessing/"
    d = joinpath("ValidationData", path)
    parsim(d,
            joinpath("ValidationRuns", path,
            @getmodel(SossModels.smring1wfVACP)),
            SossModels.smring1wfVACP
          )
    parsim(d,
          joinpath("ValidationRuns", path,
          @getmodel(SossModels.smring2wfVACP)),
          SossModels.smring2wfVACP
        )
    parsim(d,
        joinpath("ValidationRuns", path,
        @getmodel(SossModels.smring3wfVACP)),
        SossModels.smring3wfVACP
      )
    parsim(d,
      joinpath("ValidationRuns", path,
      @getmodel(SossModels.smring4wfVACP)),
      SossModels.smring4wfVACP
    )
end
