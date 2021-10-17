
function create_paths(band::String, day::String)
  ipath = joinpath("snapshot_120", "noisefrac0.02")
  pdirs = filter(x->occursin(band, x)&&occursin(day, x), readdir(datadir("sims/GRMHDCal"), join=true))
  ddirs = joinpath.(pdirs, Ref(ipath), basename.(pdirs))
  
  rdirs = projectdir.(Ref("_research"), Ref("GRMHDCal"), basename.(pdirs), Ref(ipath))
  return ddirs
end

function create_filelist(ddirs)
  return vcat(readdir.(ddirs, join=true)...)
end

function create_outdirs(filelist)
  replace.(filelist, Ref("data/sims"=>"_research")) .|> dirname .|> dirname
end
