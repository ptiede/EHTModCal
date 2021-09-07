using ROSESoss
using PyCall
using Plots
using StatsPlots
using CSV
using DataFrames
using CairoMakie
using StatsBase
include("fitscan.jl")
@pyimport corner as corn
@pyimport numpy as np



obs = ehtim.obsdata.load_uvfits("hops_3599_staticring_lo_netcal_LMTcal_JCMTcal_10s_regionIII_preimcal.uvfits")
obs_scan = obs.split_obs()
for i in 1:length(obs_scan)
    obs_scan[i].add_amp(debias=true)
    obs_scan[i].add_cphase()
end

scan=2

#########################################################3
################# Real Vis
amps = extract_amps(obs_scan[scan])
cphase = extract_cphase(obs_scan[scan])
vis = extract_vis(obs_scan[scan])

dirvis = "UniformRing/WholeData/Mring1Vis"
makedirs(dirvis)
dsum = fit_scan(dir, SossModels.mring1Vis, vis, scan, dynesty_sampler, fitgains=true)
dfvis = CSV.File(joinpath(dirvis, "Chain/nested_chain_scan-$(scan).csv"))|>DataFrame


#########################################################3
################# Real VACP
amps = extract_amps(obs_scan[scan])
cphase = extract_cphase(obs_scan[scan])
vis = extract_vis(obs_scan[scan])
dir = "UniformRing/WholeData/Mring1VACP"
makedirs(dir)
dsum = fit_scan(dir, SossModels.mring1VACP, amps, cphase, scan, dynesty_sampler, fitgains=true)
dfvacp = CSV.File(joinpath(dir, "Chain/nested_chain_scan-$(scan).csv"))|>DataFrame

########################################################################
########################################################################3
#### Simulated VACP

dirsim = "UniformRing/SynthData/Mring1VACP"
makedirs(dirsim)
#Now lets create some synthetic data
logj = ROSESoss.create_joint(SossModels.mring1VACP, amps, cphase; fitgains=true)

m = logj.model
df[end, :ma1] = 0.0
opt = copy(df[end,3:end-1])
optng = copy(df[end,11:(end-1)])
@assert opt[:ma1] == 0.0
simdata = Soss.predict(m, opt)
amps_sim = deepcopy(amps)
cp_sim = deepcopy(cphase)

#Now reset the data
amps_sim.data.amp .= simdata[:amp]
cp_sim.data.phase .= simdata[:cphase]
dsum_sim = fit_scan(dirsim, SossModels.mring1VACP, amps_sim, cp_sim, scan, dynesty_sampler, fitgains=true)
dfsim = CSV.File(joinpath(dirsim, "Chain/nested_chain_scan-$(scan).csv"))|>DataFrame


##############################3
###########################################################33
###  Visibility simulated
dirsim_vis = "UniformRing/SynthData/Mring1Vis"
makedirs(dirsim_vis)


logjvis = ROSESoss.create_joint(SossModels.mring1Vis, vis; fitgains=true)
mvis = logjvis.model
simdata_vis = Soss.predict(mvis, optng)
vis_sim = deepcopy(vis)

#Now reset the data
vis_sim.data.visr .= simdata_vis[:visr]
vis_sim.data.visi .= simdata_vis[:visi]
dsum_sim = fit_scan(dirsim_vis, SossModels.mring1Vis, vis_sim, scan, dynesty_sampler, fitgains=true)
dfsim_vis = CSV.File(joinpath(dirsim_vis, "Chain/nested_chain_scan-$(scan).csv"))|>DataFrame


##################################################################3
####################################################################3
#################3 Plots
fig = CairoMakie.density(-dfvacp[:,:mp1]*180.0/π,
                    strokecolor=(:orange),
                    strokewidth=3,
                    color=(:white, 0.0),
                    boundary=(-180.0,180.0),
                    fill=false,
                    label="Scattered Uniform Ring VACP")
CairoMakie.density!(-dfvis[:,:mp1]*180.0/π,
                    strokecolor=(:blue),
                    strokewidth=3,
                    color=(:white, 0.0),
                    boundary=(-180.0,180.0),
                    fill=false,
                    label="Scattered Uniform Ring Complex-Vis")
CairoMakie.density!(-dfsim[:,:mp1]*180.0/π,
                    strokecolor=(:green),
                    strokewidth=3,
                    color=(:white, 0.0),
                    boundary=(-180.0,180.0),
                    fill=false,
                    label="Non-scattered Uniform Ring VACP")
CairoMakie.density!(-dfsim_vis[:,:mp1]*180.0/π,
                    strokecolor=(:purple),
                    strokewidth=3,
                    color=(:white, 0.0),
                      boundary=(-180.0,180.0),
                      fill=false,
                      label="Non-scattered Uniform Ring Complex-Vis")

axislegend(position=:lt)
fig.axis.xlabel="m=1 phase (deg)"
fig
save("UniformRing/mp1_comp.png", fig)

fig = CairoMakie.hist(df[:,:ma1]*2, color=(:orange,0.5),
                      normalization=:pdf,
                      label="Scattered Uniform Ring VACP")
CairoMakie.hist!(dfsim[:,:ma1]*2, color=(:blue, 0.5),
                      normalization=:pdf,
                      label="Non-scattered Uniform Ring VACP")
CairoMakie.hist!(dfsim_vis[:,:ma1]*2, color=(:green, 0.5),
                      normalization=:pdf,
                      label="Non-scattered Uniform Ring Complex-Vis")

fig
