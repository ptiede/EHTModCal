using ROSESoss: dynesty
using Pkg; Pkg.activate(@__DIR__)

using ROSESoss
const Dists = ROSESoss.Dists
using DataFrames
using CSV
using StatsBase
using BlackBoxOptim
using CMAEvolutionStrategy
using TransformVariables
using LogDensityProblems
using ForwardDiff
using DynamicHMC
using MCMCChains
using CairoMakie
using ZigZagBoomerang
using SparseArrays
using StructArrays
include("fitscan.jl")



function makedirs(dir)
    mkpath(joinpath(@__DIR__, dir))
    mkpath(joinpath(@__DIR__, dir, "Chain"))
    mkpath(joinpath(@__DIR__, dir, "MAP"))
    mkpath(joinpath(@__DIR__, dir, "Residual"))
    mkpath(joinpath(@__DIR__, dir, "Stats"))
end

obs = ehtim.obsdata.load_uvfits("hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s_regionIII_preimcal.uvfits")
obs_scans = obs.split_obs(t_gather=60*3)
vis = extract_vis.(obs_scans)

const fwhmfac = 2*sqrt(2*log(2))


mcamp = @model u1a,v1a,u2a,v2a,u3a,v3a,u4a,v4a,errcamp,
               u1p,v1p,u2p,v2p,u3p,v3p,errcp begin
    diam ~ Dists.Uniform(25.0, 85.0)
    fwhm ~ Dists.Uniform(1.0, 40.0)
    rad = diam/2
    σ = fwhm/fwhmfac

    #First mring mode
    ma1 ~ Dists.Uniform(0.0,0.5)
    mp1 ~ Dists.Uniform(-1π,1π)
    α1 = ma1*cos(mp1)
    β1 = ma1*sin(mp1)

    #Second mring mode
    ma2 ~ Dists.Uniform(0.0,0.5)
    mp2 ~ Dists.Uniform(-1π,1π)
    α2 = ma2*cos(mp2)
    β2 = ma2*sin(mp2)

    #Fraction of floor flux
    floor ~ Dists.Uniform(0.0, 1.0)

    mring = renormed(ROSE.MRing(rad, (α1,α2), (β1,β2)), 1-floor)
    disk = renormed(stretched(ROSE.Disk(), rad, rad), floor)
    img = smoothed(mring+disk,σ)

    lcamp ~ For(eachindex(errcamp)) do i
        ca = ROSE.logclosure_amplitude(img,
                                      u1a[i], v1a[i],
                                      u2a[i], v2a[i],
                                      u3a[i], v3a[i],
                                      u4a[i], v4a[i],
                                    )
        Dists.Normal(ca, errcamp[i])
    end

    cphase ~ For(eachindex(u1cp, errcp)) do i
        mphase = ROSE.closure_phase(img,
                                    u1cp[i],
                                    v1cp[i],
                                    u2cp[i],
                                    v2cp[i],
                                    u3cp[i],
                                    v3cp[i]
                                )
        ROSE.CPVonMises(mphase, errcp[i])
    end
end

obs = ehtim.obsdata.load_uvfits("SgrAData/SgrA_ER6_ver2021-04-15/3599/hops/hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s/preimcal_pipeline_camp/snapshot60/noisefrac0.02/hops_3599_SGRA_LO_netcal_LMTcal_normalized_10s_preprocessed_snapshot60_noisefrac0.02_scan183.uvfits")
obs.add_logcamp(debias=true, count="min")
obs.add_cphase()

dcamp = ROSESoss.extract_lcamp(obs)
u1a = getdata(dcamp, :u1)
v1a = getdata(dcamp, :v1)
u2a = getdata(dcamp, :u2)
v2a = getdata(dcamp, :v2)
u3a = getdata(dcamp, :u3)
v3a = getdata(dcamp, :v3)
u4a = getdata(dcamp, :u4)
v4a = getdata(dcamp, :v4)
lcamp = getdata(dcamp, :amp)
errcamp = getdata(dcamp, :error)

dcp = ROSESoss.extract_cphase(obs)
u1cp = getdata(dcp, :u1)
v1cp = getdata(dcp, :v1)
u2cp = getdata(dcp, :u2)
v2cp = getdata(dcp, :v2)
u3cp = getdata(dcp, :u3)
v3cp = getdata(dcp, :v3)
cp = getdata(dcp, :phase)
errcp = getdata(dcp, :error)

m = mcamp(u1a=u1a,v1a=v1a,
          u2a=u2a,v2a=v2a,
          u3a=u3a,v3a=v3a,
          u4a=u4a,v4a=v4a,
          errcamp=errcamp,
          u1cp=u1cp,v1cp=v1cp,
          u2cp=u2cp,v2cp=v2cp,
          u3cp=u3cp,v3cp=v3cp,
          errcp=errcp
)

cm = m|(lcamp=lcamp, cphase=cp,)
t = xform(cm)
logdensity(cm, t(zeros(t.dimension)))





mimg = @model u, v, s1, s2, err begin
    diam ~ Dists.Uniform(25.0, 75.0)
    fwhm ~ Dists.Uniform(1.0, 40.0)
    rad = diam/2
    σ = fwhm/fwhmfac

    cij ~ Dists.Dirichlet(2, 1.0)
    ϵ ~ Dists.truncated(Dists.Normal(1.0, 0.8), 0.0, 5.0)
    size ~ Dists.truncated(Dists.Normal(50.0, 15.0), 0.0, Inf)
    ξrim ~ Dists.Uniform(-π/2, π/2)
    x0 ~ Dists.Normal(0.0, 50.0)
    y0 ~ Dists.Normal(0.0, 50.0)
    rimg = shifted(rotated(stretched(RImage(reshape(cij, 2,1), SqExpKernel(ϵ)), size, size), ξrim), x0,y0)


    #First mring mode
    ma1 ~ Dists.Uniform(0.0,0.5)
    mp1 ~ Dists.Uniform(0.0,2π)
    α1 = ma1*cos(mp1)
    β1 = ma1*sin(mp1)

    #=
    #Second mring mode
    ma2 ~ Dists.Uniform(0.0,0.5)
    mp2 ~ Dists.Uniform(0.0,2π)
    α2 = ma2*cos(mp2)
    β2 = ma2*sin(mp2)

    #Stretch
    #τ ~ Dists.Uniform(0.0, 0.99)
    #ξτ ~ Dists.Uniform(-π/2, π/2)
    #scx = sqrt(1-τ)
    #scy = 1/sqrt(1-τ)
    =#
    #Total flux
    f ~ Dists.Uniform(0.5, 5.0)

    #Fraction of floor flux
    floor ~ Dists.Uniform(0.0, 1.0)



    #mring = renormed(rotated(stretched(ROSE.MRing(rad, (α1,), (β1,)), scx, scy), ξτ), f-f*floor)
    mring = renormed(ROSE.MRing(rad, (α1,), (β1,)), f-f*floor)
    img = smoothed(mring,σ) + renormed(rimg, f*floor)
    #img = renormed(rimg, f)

    aAZ ~ Dists.truncated(Dists.LogNormal(0.0, 0.1), 0.0, 1.02)
    aJC ~ Dists.truncated(Dists.LogNormal(0.0, 0.1), 0.0, 1.02)
    aAP ~ Dists.truncated(Dists.LogNormal(0.0, 0.1), 0.0, 1.02)
    aSM ~ Dists.truncated(Dists.LogNormal(0.0, 0.1), 0.0, 1.02)
    aAA ~ Dists.truncated(Dists.LogNormal(0.0, 0.1), 0.0, 1.02)
    aLM ~ Dists.truncated(Dists.LogNormal(0.0, 0.2), 0.0, 1.02)
    aSP ~ Dists.truncated(Dists.LogNormal(0.0, 0.1), 0.0, 1.02)
    aPV ~ Dists.truncated(Dists.LogNormal(0.0, 0.1), 0.0, 1.02)
    ga = (AP=aAP, AZ=aAZ, JC=aJC, SM=aSM, AA=aAA, LM=aLM, SP=aSP, PV=aPV)

    pAP ~ Dists.VonMises(0.0, 1π)
    pAZ ~ Dists.VonMises(0.0, 1π)
    pJC ~ Dists.VonMises(0.0, 1π)
    pSM ~ Dists.VonMises(0.0, 1π)
    pAA ~ Dists.VonMises(0.0, 1π)
    pLM ~ Dists.VonMises(0.0, 1π)
    pSP ~ Dists.VonMises(0.0, 1π)
    pPV ~ Dists.VonMises(0.0, 1π)
    gp = (AP=pAP, AZ=pAZ, JC=pJC, SM=pSM, AA=pAA, LM=pLM, SP=pSP, PV=pPV)


    vis = ROSE.visibility.(Ref(img), u,v)
    visr ~ For(eachindex(u, v)) do i
        Δθ = gp[s1[i]] - gp[s2[i]]
        s,c = sincos(Δθ)
        g1 = ga[s1[i]]
        g2 = ga[s2[i]]
        mamp = g1*g2*(real(vis[i])*c + imag(vis[i])*s)
        Dists.Normal(mamp, err[i])
    end

    visi ~ For(eachindex(u, v)) do i
        Δθ = gp[s1[i]] - gp[s2[i]]
        s,c = sincos(Δθ)
        g1 = ga[s1[i]]
        g2 = ga[s2[i]]
        mamp = g1*g2*(-real(vis[i])*s + imag(vis[i])*c)
        Dists.Normal(mamp, err[i])
    end
end

Soss.xform(d::Dists.Dirichlet, _data::NamedTuple) = UnitSimplex(length(d.alpha))


lj = create_joint(mimg, vis[3]; fitgains=true)


m = lj.model
cm = m|(lj.data)
t = lj.transform

ℓ(x) = logdensity(cm, TransformVariables.transform(t, x))
P = LogDensityProblems.TransformedLogDensity(t, ℓ)
∇P = LogDensityProblems.ADgradient(:ForwardDiff, P)

ℓt(x) = LogDensityProblems.logdensity(P, x)

function construct_ranges(t)
    ranges = Tuple{Float64,Float64}[]
    for i in 1:length(t)
        if typeof(t[i]) <: TransformVariables.ScaledShiftedLogistic
            push!(ranges, (-5.0,5.0))
        elseif typeof(t[i]) <: TransformVariables.ShiftedExp
            push!(ranges, (0.0, 10.0))
        elseif typeof(t[i]) <: TransformVariables.UnitSimplex
            tmp = [(-5.0,5.0) for i in 1:(t[i].n-1)]
            push!(ranges, tmp...)
        elseif typeof(t[i]) <: TransformVariables.Identity
            push!(ranges, (-10.0, 10.0))
        else
            throw("Tranform $(typeof(t[i])) not found")
        end
    end
    return ranges
end
sranges = construct_ranges(t.transformations)



function threaded_opt(nopt, f, ranges)
    pos = Vector{Vector{Float64}}(undef, nopt)
    minf = Vector{Float64}(undef, nopt)
    lower = first.(ranges)
    upper = last.(ranges)
    Threads.@threads for i in 1:nopt
        start = lower .+ rand(length(ranges)).*(upper .- lower)
        res = minimize(f, start, 1., lower=lower, upper=upper)
        pos[i] = xbest(res)
        minf[i] = fbest(res)
    end
    idxs = sortperm(minf)
    return pos[idxs], minf[idxs]
end
pos, minf = threaded_opt(72, x->-ℓ(x), sranges)

maxp = TransformVariables.transform(t, pos[5])
ff = TransformVariables.transform(t, pos[2])
img = Soss.predict(m, maxp)[:img]
sim = intensitymap(img, 128, 128, 120.0,120.0)
image(imagepixels(sim)..., sim', axis=(aspect=1,xreversed=true), colormap=:viridis)



reporter = DynamicHMC.LogProgressReport(step_interval=10)
results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, 5000;
                           warmup_stages = default_warmup_stages(doubling_stages=8),
                           initialization=(q=pos[1],), reporter = reporter)
chain = TransformVariables.transform.(Ref(t), first(results))

imgs = getproperty.(Soss.predict.(Ref(m), chain[1:10:end]), Ref(:img))
sims = intensitymap.(imgs, 128, 128, 200.0, 200.0)
image(imagepixels(sims[1])..., sims[490]', axis=(aspect=1, xreversed=true), colormap=:viridis)


function make()
    r0 = 9.5e-11*1e6*3600*180/π
    floor=0.01
    ma1 = 0.17
    mp1 = -2.27
    b1,a1 = ma1.*sincos(mp1)
    ma2 = 0.25
    mp2 = -1.82
    b2,a2 = ma2.*sincos(mp2)
    ma3 = 0.47
    mp3 = -1.17
    b3,a3 = ma3.*sincos(mp3)
    ma4 = 0.28
    mp4 = 0.0273
    b4,a4 = ma4.*sincos(mp4)
    τ = 0.46
    scx = 1/sqrt(1-τ)
    scy = sqrt(1-τ)
    ξτ = π-0.77
    σ = 5.8e-11*180.0/π*3600*1e6
    m = @chain begin
        MRing(r0, (a1,a2,a3,a4),(b1,b2,b3,b4))
        stretched(scx,scy)
        rotated(ξτ)
        smoothed(σ)
    end
    return m
end
