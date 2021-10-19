using Pkg
Pkg.activate(@__DIR__)

using RobustAdaptiveMetropolisSampler
#Do I really need to do this? I guess so because of SO that should be fixed in 1.7
sample(x->-sum(abs2, x), RAM(zeros(10), 1.0), 100000, show_progress=true)
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


function quickprocess(mdir; tmin=0.0, tmax=24.0, nbatch=1_000, ckpt_stride=50_000, nsteps=2_000_000)
    mins, maxs, wrapped, quants, labels = parsechainpath(mdir)
    outdir = joinpath(mdir, "ChainHA_2")

    process(joinpath(mdir, "chain.h5"), mins, maxs, wrapped, quants, labels, outdir; nbatch, ckpt_stride, tmin, tmax, nsteps)

end

function create_lklhd(cfile, mins, maxs, wrapped, quants; tmin=0.0, tmax=24.0, nbatch=1000)
    chain = ChainH5(cfile, quants)
    chainall = restricttime(chain, tmin, tmax)
    l = BatchStackerLklhd(chainall, mins, maxs, wrapped, nbatch)
    prior = (μ = Product(Uniform.(mins, maxs)), σ = Product(Uniform.(0.01, maxs .- mins)))

    return l, prior, keys(chainall)
end

function write_results(file,tv, logp, keys; stride=1)
    echain = TupleVector(tv[1:stride:end]) #clip the first 50% of the chain
    df = _mkdf(echain, keys)
    df.logp = logp[1:stride:end]
    df|>CSV.write(file)
end


function process(cfile, mins, maxs, wrapped, quants, labels, outdir; ckpt_stride=50_000, tmin=0, tmax=24.0, stride=10, nbatch=1_000, nsteps=2_000_000)
    #diameter
    mkpath(outdir)
    l, prior, keys = create_lklhd(cfile, mins, maxs, wrapped, quants; tmin=tmin, tmax=tmax, nbatch=nbatch)
    p0, fmap = optimize(l, prior, ROSESoss.BBO(;maxevals=25_000))
    if ckpt_stride > nsteps
        println("Checkpoint stide > nsteps resetting")
        ckpt_stride = nsteps
    end
    Minit = 0.01#[0.01*(maxs .- mins)..., 0.001*(maxs .- mins)...]
    smplr = RAM(p0, Minit)
    nbatches = nsteps÷ckpt_stride
    res = @timed sample(l, prior, smplr, ckpt_stride; show_progress=false, output_log_probability_x=true)
    tv = res.value[1]
    state = res.value[3]
    logp = res.value[4]
    println("Done batch 1 this took $(res.time) seconds")
    println("I am guessing this entire run will take $(res.time*nbatches/3600.0) hours to finish")
    write_results(joinpath(outdir, replace(basename(cfile), ".h5"=>"_ha.csv")),tv, logp, keys)
    for i in 2:nbatches
        println("On batch $i/$nbatches")
        tv_b,stats_b,state_b,logp_b = sample(l, prior, state, ckpt_stride; show_progress=false, output_log_probability_x=true)
        #extend results and set up next run
        push!(tv, tv_b...)
        push!(logp, logp_b...)
        state = state_b
        write_results(joinpath(outdir, replace(basename(cfile), ".h5"=>"_ha.csv")), tv, logp, keys, stride=10)
    end

    return nothing
end
