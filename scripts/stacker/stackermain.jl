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
using Serialization

include(joinpath(@__DIR__,"build_prob.jl"))
include(joinpath(@__DIR__,"converters.jl"))
include(joinpath(@__DIR__, "inference.jl"))


function quickprocess(mdir; restart=false, tmin=0.0, tmax=24.0, nbatch=1_000, ckpt_stride=100_000, nsteps=2_000_000)
    mins, maxs, wrapped, quants, labels = parsechainpath(mdir)
    outdir = joinpath(mdir, "ChainHA_3")

    process(joinpath(mdir, "chain.h5"), mins, maxs, wrapped, quants, labels, outdir; restart, nbatch, ckpt_stride, tmin, tmax, nsteps)

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

function readckpt(ckptfile)
    return deserialize(ckptfile)
end

function saveckpt(ckptfile, smplr, tv, logp)
    return serialize(ckptfile, (;smplr, tv, logp))
end


function process(cfile, mins, maxs, wrapped, quants, labels, outdir; restart=false, ckpt_stride=50_000, tmin=0, tmax=24.0, stride=10, nbatch=1_000, nsteps=2_000_000)
    #diameter

    #Create output stuff
    mkpath(outdir)
    ckptfile = joinpath(outdir, replace(basename(cfile), ".h5"=>"_ckpt.jls"))
    chainfile = joinpath(outdir, replace(basename(cfile), ".h5"=>"_ha_trunc.csv"))
    nbatches = nsteps÷ckpt_stride

    l, prior, keys = create_lklhd(cfile, mins, maxs, wrapped, quants; tmin=tmin, tmax=tmax, nbatch=nbatch)
    if !restart || !isfile(ckptfile)

        p01, fmap1 = optimize(l, prior, ROSESoss.BBO(;maxevals=40_000))
        p02, fmap2 = optimize(l, prior, ROSESoss.BBO(;maxevals=40_000))
        println("After 2 runs the estiamte maps are $(fmap1) and $(fmap2)")
        #Run twice and take better
        if fmap1 > fmap2
            p0 = p01
        else
            p0 = p02
        end
        if ckpt_stride > nsteps
            println("Checkpoint stide > nsteps resetting")
            ckpt_stride = nsteps
        end
        Minit = 0.001#[0.01*(maxs .- mins)..., 0.001*(maxs .- mins)...]
        smplr = RAM(p0, Minit)
        res = @timed sample(l, prior, smplr, ckpt_stride; show_progress=false, output_log_probability_x=true)
        tv = res.value[1]
        state = res.value[3]
        logp = res.value[4]
        # Starting checkpoint
        saveckpt(ckptfile, state, tv, logp)
        println("Done batch 1 this took $(res.time) seconds")
        println("I am guessing this entire run will take $(res.time*nbatches/3600.0) hours to finish")
        write_results(chainfile,tv, logp, keys)
        nstart=2
    else
        state, tv, logp = readckpt(ckptfile)
        nstart = length(tv)÷ckpt_stride + 1
        @info "Reading in checkpoint file $(ckptfile)"
        @info "According to checkpoint I am on batch $(nstart)"
    end

    for i in nstart:nbatches
        println("On batch $i/$nbatches")
        # We are early so lets speed up the learning
        if i < nbatches÷2
            state = RAM(state.x, state.M.mat; step=1)
        end
        tv_b,state_b,state_b,logp_b = sample(l, prior, state, ckpt_stride; show_progress=false, output_log_probability_x=true)

        #extend results and set up next run
        push!(tv, tv_b...)
        push!(logp, logp_b...)
        state = state_b
        println("Writing checkpoint")
        saveckpt(ckptfile, state, tv, logp)
        write_results(chainfile, tv, logp, keys, stride=stride)
    end

    return nothing
end
