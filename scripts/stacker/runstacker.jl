using DrWatson
@quickactivate
using EHTModelStacker
using CairoMakie
using PyCall
using TupleVectors
using HDF5
@pyimport dynesty

using ParameterHandling
using StatsBase
using CSV, DataFrames
using Distributions
using KernelDensity
using HypercubeTransform
using StaticArrays
using Bijectors
using RobustAdaptiveMetropolisSampler
using EHTModelStacker: MvNormal2D, MvNormalFast
using ROSESoss
using ArraysOfArrays
using NamedTupleTools

include(joinpath(@__DIR__,"build_prob.jl"))
include(joinpath(@__DIR__,"converters.jl"))
include(joinpath(@__DIR__,"plots.jl"))


function average_chain_diag(build_model, mins, maxs, wrapped, chain; nlive=1500)
    @assert length(mins)==length(maxs)
    @assert length(wrapped)==length(maxs)
    priors = (μ = Product(Uniform.(mins, maxs)), σ = Product(Uniform.(0.0, maxs .- mins)))
    t = ascube(priors)
    p0 = HypercubeTransform.transform(t, fill(0.5, dimension(t)))
    fp0, unflatten = ParameterHandling.flatten(p0)
    lp = build_loglklhd(build_model, chain, mins, maxs, wrapped)∘unflatten
    pt(x) = first(ParameterHandling.flatten(HypercubeTransform.transform(t, x)))

    # return lp, pt, unflatten
    sampler = dynesty.NestedSampler(lp, pt, dimension(t), nlive=nlive)
    sampler.run_nested()
    res = sampler.results
    samples, weights = res["samples"], exp.(res["logwt"] .- res["logz"][end])
    tv = TupleVector([unflatten(Vector(x)) for x in eachrow(samples)])
    echain = TupleVector(sample(tv, Weights(weights), 4000))
    return echain

end


function quickprocess(cfile; plotres=true, tmin=0, tmax=1e30, nlive=1500)
    mins, maxs, wrapped, quants, labels = parsechainpath(cfile)
    outdir = joinpath(cfile, "ChainHA_W")

    process(joinpath(cfile, "chain.h5"), mins, maxs, wrapped, quants, labels, outdir; plotres, tmin, tmax, nlive)

end




function process(cfile, mins, maxs, wrapped, quants, labels, outdir; plotres=true, tmin=0, tmax=1e30, nlive=1500)
    #diameter
    mkpath(outdir)
    chain = ChainH5(cfile, quants, 100)
    chainall = restricttime(chain, tmin, tmax)
    echain = average_chain_diag(build_model_all, mins, maxs, wrapped, chainall; nlive)
    df = _mkdf(echain, keys(chainall))
    df|>CSV.write(joinpath(outdir, replace(basename(cfile), ".h5"=>"_ha.csv")))
    if plotres
        for (i,q) in enumerate(quants)
            c = getparam(chainall, q)
            ec = TupleVector((μ=getindex.(echain.μ,i), σ=getindex.(echain.σ,i)))
            fig = plot_stacker_summary(c, ec, labels[i], mins[i], maxs[i])
            save(joinpath(outdir, replace(basename(cfile), ".h5"=>"_$q.png")), fig)
        end
    end
    return nothing
end
