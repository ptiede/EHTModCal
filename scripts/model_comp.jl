using ROSE: visibility_amplitude, uvtriangle
using ROSE
using CairoMakie
import ROSESoss

using Soss
using PyCall
@pyimport ehtim

obs = ehtim.obsdata.load_uvfits("hops_3599_SGRA_HI_netcal_10s.uvfits")
obsavg = obs.avg_coherent(12.0*60.0)
obsavg.add_amp(debias=true)
obsavg.add_cphase()
amp = ROSESoss.extract_amps(obsavg)
cphase = ROSESoss.extract_cphase(obsavg)

uvtriangle = ROSE.uvtriangle.(cphase.data)

u = getdata(amp, :u)
v = getdata(amp, :v)

fig = CairoMakie.scatter(u, v)
CairoMakie.scatter!(-u, -v)

m = smoothed(renormed(ROSE.ConcordanceCrescent(65.0/2, 30.0/2, 0.0, 0.5), 1.0),0.5)

sim = intensitymap(m, 128, 128, 200.0,200.0)
image(imagepixels(sim)..., sim', colormap=:viridis, axis=(aspect=1, xreversed=true))

err = getdata(amp, :error)
ma1 = 0.25
mp1 = -π/2
b1, a1 = ma1.*sincos(mp1)
mr = smoothed(MRing(22.5, (a1,), (b1,)), 5.0)
simr = intensitymap(mr, 256, 256, 200.0,200.0)
image(imagepixels(simr)..., simr', colormap=:viridis, axis=(aspect=1, xreversed=true))

x = -120.0:0.1:120.0
fl = lines(imagepixels(simr)[1], simr[128,:], color=:red)
lines!(imagepixels(sim)[1], sim[64,:], color=:blue)
fl

ampfake = ROSE.visibility_amplitude.(Ref(m), u, v)
cpfake = [ROSE.closure_phase(m, uvtri...) for uvtri in uvtriangle]

logj = ROSESoss.create_joint(ROSESoss.SossModels.smring2wfVACP, amp, cphase, fitgains=false)

sossm = logj.model

csossm = sossm | (amp=ampfake, cphase=cpfake, AP = 1.0, AZ = 1.0, JC = 1.0, SM = 1.0, AA = 1.0, LM = 1.0, SP = 1.0, PV = 1.0)

res = ROSESoss.dynesty_sampler(csossm)
chain = first(res);
echain = sample(chain, Weights(vec(chain[:weights])), 3000)
lines(vec(echain[:mp1]))

dfchain = DataFrame(echain)
r1 = NamedTuple(dfchain[end, 3:end-1])

mfit = Soss.predict(csossm, r1)[:img]
simfit = intensitymap(mfit, 256, 256, 200.0,200.0)

fl = lines(imagepixels(sim)[1], sim[64,:], color=:blue)
lines!(imagepixels(simfit)[1], simfit[128,:])
fl

fig = CairoMakie.scatter(hypot.(u,v), ampfake, color=:blue)
CairoMakie.scatter!(hypot.(u,v), ROSE.visibility_amplitude.(Ref(mfit), u, v), color=:red)
fig