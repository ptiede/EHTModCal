using DrWatson
@quickactivate
using ROSESoss
using TupleVectors
using HypercubeTransform
using StatsBase

obs = ehtim.obsdata.load_uvfits(datadir("sims/ValidationData/ehtim_geometric_model/obs_mring_synthdata_thermal_only.uvfits"))
obs_avg = obs.copy().avg_coherent(180.0)
obs_avg.add_logcamp(debias=true, count="min")
obs_avg.add_cphase(count="min-cut0bl")


dlca = ROSESoss.extract_lcamp(obs_avg)
dcp = ROSESoss.extract_cphase(obs_avg)

cm = create_joint(ROSESoss.mringwfloor(N=3,), dlca, dcp)
tc = ascube(cm)

tv, stats = sample(DynestyStatic(nlive=1000), cm)

tve = sample(tv, Weights(tv.weights), 5000) |> TupleVector