using Printf

function create_p6dirs(dataset, time)
  prepath = "/gpfs/ptiede/EHTModCal/_research/P6ModCal/"
  dstr = @sprintf "dataset%03d" dataset
  if dataset < 36
    path = joinpath(prepath, "uvfits_MoD_precal_stage1")
  elseif dataset < 119
    path = joinpath(prepath, "uvfits_MoD_precal_stage2")
  end
  return readdir(joinpath(path,dstr, time, "noisefrac0.02"), join=true)
end
