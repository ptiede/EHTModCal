using DrWatson
@quickactivate

using ROSESoss
import MeasureTheory as MT
import Distributions as Dists



obs = ehtim.obsdata.load_uvfits("/home/ptiede/Research/Projects/RefractiveScatteringThemis/TestDataGen/ER5/lo/hops_3601_M87+netcal.uvfits")
obs.add_scans()
obsavg = obs.avg_coherent(0.0, scan_avg=true)
obsavg.add_logcamp(debias=true, count="min")
obsavg.add_cphase(count="min-cut0bl")
obsavg.add_amp(debias=true)


dlca = ROSESoss.extract_lcamp(obsavg)
dcp = ROSESoss.extract_cphase(obsavg)
dva = ROSESoss.extract_amps(obsavg)

cm = create_joint(ROSESoss.mring(N=1,), dlca, dcp)
dy = DynestyStatic(nlive=200)
tvd, statsd = sample(dy, cm; print_progress=true)


cmva = create_joint(ROSESoss.mring(N=1,), dva, dcp)
tvd, statsd, sample(dy, cmva; print_progress=true)
