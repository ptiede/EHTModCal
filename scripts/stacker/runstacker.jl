using DrWatson
@quickactivate
using EHTModelStacker
using CairoMakie
using PyCall
using TupleVectors
using HDF5
using ROSESoss

using ParameterHandling
using StatsBase
using CSV, DataFrames
using Distributions
using KernelDensity
using HypercubeTransform
using StaticArrays
using RobustAdaptiveMetropolisSampler
using EHTModelStacker: MvNormal2D, MvNormalFast, MvUniform
using ArraysOfArrays
using NamedTupleTools
using BlackBoxOptim

include(joinpath(@__DIR__,"build_prob.jl"))
include(joinpath(@__DIR__,"converters.jl"))
include(joinpath(@__DIR__,"plots.jl"))
include(joinpath(@__DIR__, "inference.jl"))





function quickprocess(mdir; plotres=true, tmin=0, tmax=1e30, nsteps=500_000)
    mins, maxs, wrapped, quants, labels = parsechainpath(mdir)
    outdir = joinpath(mdir, "ChainHA_W")

    process(joinpath(mdir, "chain.h5"), mins, maxs, wrapped, quants, labels, outdir; plotres, tmin, tmax, nsteps)

end




function process(cfile, mins, maxs, wrapped, quants, labels, outdir; plotres=true, tmin=0, tmax=24.0, nbatch=2_00, nsteps=2_000_000)
    #diameter
    mkpath(outdir)
    chain = ChainH5(cfile, quants)
    chainall = restricttime(chain, tmin, tmax)
    l = BatchStackerLklhd(chain, mins, maxs, wrapped, nbatch)
    prior = (μ = Product(Uniform.(mins, maxs)), σ = Product(Uniform.(0.01, maxs .- mins)))
    xbest, lmap = ROSESoss.optimize(l, prior, ROSESoss.BBO(;maxevals=10_000, tracemode=:silent))
    println("lmap is $lmap")
    Minit = [0.01*(maxs .- mins)..., 0.001*(maxs .- mins)...]
    tv, stats, M0 = sample(l, prior, RAM(M0=Minit, n=nsteps, show_progress=false), xbest)
    echain = TupleVector(tv[end÷2:end]) #clip the first 50% of the chain
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
