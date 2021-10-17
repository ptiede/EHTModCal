using Pkg
Pkg.activate(@__DIR__)

using RobustAdaptiveMetropolisSampler
#Do I really need to do this? I guess so because of SO that should be fixed in 1.7
RAM_sample(x->-sum(abs2, x), zeros(10), 1.0, 100000, show_progress=true)
using ROSESoss
using BlackBoxOptim
using ParameterHandling
using Distributions
using DataFrames
using CSV
using HypercubeTransform
using StatsBase
using EHTModelStacker
using TupleVectors

include(joinpath(@__DIR__,"build_prob.jl"))
include(joinpath(@__DIR__,"converters.jl"))
include(joinpath(@__DIR__, "inference.jl"))


function quickprocess(mdir; tmin=0.0, tmax=24.0, nbatch=1_000, nsteps=2_000_000)
    mins, maxs, wrapped, quants, labels = parsechainpath(mdir)
    outdir = joinpath(mdir, "ChainHA_2")

    process(joinpath(mdir, "chain.h5"), mins, maxs, wrapped, quants, labels, outdir; nbatch,  tmin, tmax, nsteps)

end

function create_lklhd(cfile, mins, maxs, wrapped, quants; tmin=0.0, tmax=24.0, nbatch=1000)
    chain = ChainH5(cfile, quants)
    chainall = restricttime(chain, tmin, tmax)
    l = BatchStackerLklhd(chainall, mins, maxs, wrapped, nbatch)
    prior = (μ = Product(Uniform.(mins, maxs)), σ = Product(Uniform.(0.01, maxs .- mins)))

    return l, prior, keys(chainall)
end

function write_results(file,tv, keys)
    echain = TupleVector(tv[end÷2:end]) #clip the first 50% of the chain
    df = _mkdf(echain, keys)
    df|>CSV.write(file)
end


function process(cfile, mins, maxs, wrapped, quants, labels, outdir; tmin=0, tmax=24.0, nbatch=1_000, nsteps=2_000_000)
    #diameter
    mkpath(outdir)
    l, prior, keys = create_lklhd(cfile, mins, maxs, wrapped, quants; tmin=tmin, tmax=tmax, nbatch=nbatch)
    p0, fmap = optimize(l, prior, ROSESoss.BBO(;maxevals=5_000))
    Minit = 0.01#[0.01*(maxs .- mins)..., 0.001*(maxs .- mins)...]
    tv,stats,M0 = StatsBase.sample(l, prior, RAM(M0=Minit, n=nsteps, show_progress=true), p0)
    write_results(joinpath(outdir, replace(basename(cfile), ".h5"=>"_ha.csv")),tv, keys)
    return nothing
end
