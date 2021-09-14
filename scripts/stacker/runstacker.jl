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


function build_model_LN(θ::NamedTuple, min, max)
    μ, σ = θ.μ, θ.σ
    println(θ.σ)
    transition = LogNormal(log(μ), log(σ))
    prior = Uniform(min, max)
    return SnapshotWeights(transition, prior)
end

function build_model_N(θ::NamedTuple, min, max)
    μ, σ = θ.μ, θ.σ
    transition = Normal(μ, σ)
    prior = Uniform(min, max)
    return SnapshotWeights(transition, prior)
end

function build_model_all(θ::NamedTuple, mins, maxs)
    μ, Σ = θ.μ, θ.σ
    transition = MvNormal(μ, Σ)
    prior = product_distribution([Uniform(mins[i], maxs[i]) for i in eachindex(mins, maxs)])
    return SnapshotWeights(transition, prior)
end


function build_loglklhd(build_model, chain, min, max)
    return function (θ::NamedTuple)
        ws = build_model(θ, min, max)
        return lpdf(ws, chain)
    end
end

function build_prior_transform(minp, maxp)
    function (p)
        pμ = minp + (maxp-minp)*p[1]
        pσ = (maxp-minp)*p[2]
        return (pμ, pσ)
    end
end

function movie_uncert(chainfiles, mfile, outname)
    mvals,times = h5open(chainfiles, "r") do fid
        times = read(fid["time"])
        tt = String.(keys(fid["params"]["scan1"]))
        pnames = filter(x->!((x=="weights")||(x=="chain")||(x=="iteration")||(x=="diamdb")), tt)
        mvals = NamedTuple{Tuple(Symbol.(pnames))}[]
        for i in 1:length(times)
            tmp = Float64[]
            dp = read(fid["params"]["scan$i"])
            for k in pnames
                push!(tmp, mean(dp[k]))
            end
            push!(mvals, NamedTuple{Tuple(Symbol.(pnames))}(Tuple(tmp)))
        end
        return mvals,times
    end
    model = Soss.prior(last(parsemodel(mfile)), :amp, :cphase)
    m = model()
    fig = Figure(resolution=(400,400))
    ax = Axis(fig[1,1], xlabel="RA (μas)", ylabel="DEC (μas)", xreversed=true)
    record(fig, outname, 1:length(mvals); framerate=5) do i
        p = mvals[i]
        mim = Soss.predict(m, p)[:img]
        img = intensitymap(mim, 128,128,160.0,160.0)
        ax.title="Time = $(times[i])"
        image!(ax, imagepixels(img)..., img', colormap=:afmhot)
    end
end

function frames_uncert(chainfiles, mfile, time, nsamples, outname)
    mvals,time = h5open(chainfiles, "r") do fid
        times = read(fid["time"])
        itime = findfirst(t->t>time, times)
        tt = String.(keys(fid["params"]["scan$itime"]))
        dp = read(fid["params"]["scan$itime"])
        pnames = filter(x->!((x=="weights")||(x=="chain")||(x=="iteration")||(x=="diamdb")), tt)
        mvals = NamedTuple{Tuple(Symbol.(pnames))}[]
        for i in 1:nsamples
            tmp = Float64[]
            for k in pnames
                push!(tmp, dp[k][i])
            end
            push!(mvals, NamedTuple{Tuple(Symbol.(pnames))}(Tuple(tmp)))
        end
        return mvals,time
    end
    model = Soss.prior(last(parsemodel(mfile)), :amp, :cphase)
    m = model()
    fig = Figure(resolution=(400,400))
    ax = Axis(fig[1,1], xlabel="RA (μas)", ylabel="DEC (μas)", xreversed=true)
    record(fig, outname, 1:length(mvals); framerate=5) do i
        p = mvals[i]
        mim = Soss.predict(m, p)[:img]
        img = intensitymap(mim, 128,128,160.0,160.0)
        ax.title="Time = $(time)"
        image!(ax, imagepixels(img)..., img', colormap=:afmhot)
    end
end

function average_chain(build_model, minp, maxp, chain)
    p0 = (μ=(maxp-minp/2)+minp, σ=(maxp-minp)/2)
    fp0, unflatten = ParameterHandling.flatten(p0)
    lp = build_loglklhd(build_model, chain, minp, maxp)∘unflatten
    pt = build_prior_transform(minp, maxp)

    fu = Dict("min_ncall"=>100, "min_eff"=>0.2)
    sampler = dynesty.NestedSampler(lp, pt, 2, first_update=fu)
    sampler.run_nested()
    res = sampler.results
    samples, weights = res["samples"], exp.(res["logwt"] .- res["logz"][end])
    tv = TupleVector((μ=samples[:,1], σ=samples[:,2]))
    echain = TupleVector(sample(tv, Weights(weights), 2000))
    return echain
end

function process(cfile, minp, maxp, quant, label, outdir; plotres=true)
    #diameter
    mkpath(outdir)
    chain = ChainH5(cfile, quant)
    echain = average_chain(build_model_N, minp, maxp, chain)
    echain|>DataFrame|>CSV.write(joinpath(outdir, replace(basename(cfile), ".h5"=>"_ha_$quant.csv")))
    if plotres
      fig = Figure(resolution=(600, 800))
      ax = Axis(fig, ylabel="Probability Density")
      axsa = Axis(fig)
      kμ = kde(echain.μ)
      lower = fill(0.0, length(kμ.x))
      upper = kμ.density
      lines!(ax, kμ.x, upper, color=(:black, 0.2))
      band!(ax, kμ.x, lower, upper, color=(:blue, 0.5), label="μ $label")
      ridge(axsa, chain)
      fig[1,1] = ax
      fig[2,1] = axsa
      rowsize!(fig.layout, 2, Relative(3/4))
      kσ = kde(echain.σ)
      errorbars!(ax, mean(echain.μ), fill(mean(kμ.density), 1), 2.0*quantile(echain.σ, 0.5), whiskerwidth = 10,
              direction = :x, label="σ $label")

      linkxaxes!(axsa, ax)
      hidexdecorations!(ax, grid=false, ticks=false)
      axsa.xlabel = label
      axislegend(ax, framevisible=false)
      xlims!(ax, minp, maxp)
      save(joinpath(outdir, replace(basename(cfile), ".h5"=>"_$quant.png")), fig)
    end
    return summarize(echain)
end

function allprocess(cfile, outdir; plotres=true)
    #diameter
    mkpath(outdir)
    sdiam = process(cfile, 0.0, 100.0, "diam", "diameter (μas)", outdir; plotres=plotres)
    #width
    swidth = process(cfile, 0.0, 40.0, "fwhm", "width (μas)", outdir; plotres=plotres)
    #mp1
    spa = process(cfile, -π, π, "mp1", "PA (rad)", outdir; plotres=plotres)
    #floor
    sfloor = process(cfile, 0.0, 1.0, "floor", "floor flux", outdir; plotres=plotres)
    return (diam_μ = sdiam.μ, diam_σ = sdiam.σ,
            width_μ = swidth.μ, width_σ = swidth.σ,
            pa_μ = spa.μ, pa_σ = spa.σ,
            floor_μ = sfloor.μ, floor_σ = sfloor.σ,
           )
end

function main(dir)
    cfiles = filter(endswith(".h5"), readdir(dir))
    for c in cfiles
        name = first(splitext(c))
        println("On chain $c")
        outdir = mkpath(joinpath(dir, name*"_ha"))
        allprocess(joinpath(dir,c), outdir)
    end
end

function ridge(ax, chain::ChainH5; mod=(x->x), kdekwargs=(),)
    times = chain.times
    height = (times[end]-times[1])/length(times)*3.0
    my = map(mod, chain.chain)
    cls = to_colormap(:viridis, length(times))
    for i in 1:length(times)
        k = kde(my[i]; kdekwargs...)
        lower = fill(times[i], length(k.x))
        upper = times[i] .+ k.density*height/maximum(k.density)
        lines!(k.x, upper, color=(:black, 0.2))
        band!(k.x, lower, upper, color=(cls[i], 0.5))
    end
    #band!(ax, -180:1:180, 13.4,13.64, color=(:grey, 0.4))
    ax.ylabel = "Time UTC"
    #ax.yticks = (range(times[1],times[end-1], length=10))
    return ax
end

function plotridgemp1(cfile)
    chain = ChainH5(cfile, "mp1")
    fig = Figure(resolution=(400,800))
    ax = Axis(fig[1,1], xlabel="Position Angle  (deg)")
    ridge(ax, chain, mod=(x->-180.0/π*x), kdekwargs=(boundary=(-180,180),))
    #CairoMakie.ylims!(12.5,14.2)
    return fig
end


function make_plot(day)
    dirs = filter(isdir, readdir(joinpath("RunsCal/SgrARuns/BestTimes",
                                day,
                                "preimcal_camp/snapshot60"), join=true))
    cfiles =  joinpath.(dirs, Ref("chain.h5"))
    fig = Figure(resolution=(200*length(cfiles),600))
    axes = Axis[]
    for (i,c) in pairs(cfiles)
        ax = Axis(fig[1,i],
                  xlabel="Position Angle  (deg)",
                  title=first(split(splitpath(c)[end-1],"wf"))
            )
        chain = ChainH5(c, "mp1")
        ridge(ax, chain, mod=(x->-180.0/π*x), kdekwargs=(boundary=(-180,180),))
        push!(axes,ax)
        if i >1
            hideydecorations!(ax, grid=false, ticks=false)
        end
    end
    linkyaxes!(axes...)
    return fig,ax
end

function plotchi2(day)
    dirs = filter(isdir, readdir(joinpath("RunsCal/SgrARuns/BestTimes",
    day,
    "preimcal_camp/snapshot60"), join=true))
    cfiles =  joinpath.(dirs, Ref("merged_stats.csv"))
    dfsums = []
    sims = []
    for c in cfiles
        if occursin("smring", c)
            dfsum = CSV.File(c) |> DataFrame
            push!(dfsums, dfsum)
            push!(sims, first(split(splitpath(c)[end-1],"wf")))
        end
    end

    fig = Figure(reoslution=(300,400))
    ax1 = Axis(fig[1,1], xlabel="Time (UTC)", ylabel="χ² amp.")
    ax2 = Axis(fig[2,1], xlabel="Time (UTC)", ylabel="χ² cp")
    for (i,df) in enumerate(dfsums)
        lines!(ax1, df[!,:time], df[!,:mampchi2], label=sims[i])
        lines!(ax2, df[!,:time], df[!,:mcpchi2], label=sims[i])
    end
    linkxaxes!(ax1,ax2)
    hidexdecorations!(ax1, grid=false, ticks=false)
    leg = fig[:,2] = Legend(fig, ax1)

    return fig
end

function plotlogz(day)
    dirs = filter(isdir, readdir(joinpath("RunsCal/SgrARuns/BestTimes",
    day,
    "preimcal_camp/snapshot60"), join=true))
    cfiles =  joinpath.(dirs, Ref("chain.h5"))
    logz = Float64[]
    chains = []
    sims = []
    for (i,c) in pairs(cfiles)
        chain = ChainH5(c, "diam")
        push!(logz, sum(chain.logz))
        push!(chains, chain)
        push!(sims, first(split(splitpath(c)[end-1],"wf")))
    end
    fig = Figure(resolution=(400,400))
    ax = Axis(fig[1,1],
            xticks = (collect(1:length(logz)), sims),
            xticklabelrotation=π/4,
            ylabel="Δlogz",
            title=day)

    scatter!(ax,1:length(logz), logz .- maximum(logz))
    return fig, ax, logz

end


function comboavgplot(fchainsa, fcsa, outname)
    chainsa = ChainH5(fchainsa, "img_diam")
    csa = CSV.File(fcsa) |> DataFrame

    fig = Figure(resolution=(600, 800))
    ax = Axis(fig, ylabel="Probability Density")
    axsa = Axis(fig)
    kμ = kde(csa[:,:μ])
    lower = fill(0.0, length(kμ.x))
    upper = kμ.density
    lines!(ax, kμ.x, upper, color=(:black, 0.2))
    band!(ax, kμ.x, lower, upper, color=(:blue, 0.5), label="μ diameter")
    ridge(axsa, chainsa)
    fig[1,1] = ax
    fig[2,1] = axsa
    rowsize!(fig.layout, 2, Relative(3/4))
    kσ = kde(csa[:,:σ])
    errorbars!(ax, mean(csa[:,:μ]), fill(mean(kμ.density), 1), 2.0*quantile(csa[:,:σ], 0.5), whiskerwidth = 10,
              direction = :x, label="σ diameter")

    linkxaxes!(axsa, ax)
    hidexdecorations!(ax, grid=false, ticks=false)
    axsa.xlabel = "diameter (μas)"
    axislegend(ax, framevisible=false)
    xlims!(ax, 25.0,85.0)
    save(outname, fig)

    display(fig)
    return fig
end

function violin_all(dir, quant; mod=(x->1*x))
    dirs = filter(isdir, readdir(dir, join=true))
    files = joinpath.(dirs, basename.(dirs).*"_$quant.csv")
    fig = Figure(resolution=(1000, 600))
    ax = Axis(fig[1,1],
              xticks = (collect(1:length(files)), ["dataset $(i-1)" for i in 1:length(files)]),
              xticklabelrotation=π/2)
    for i in 1:length(files)
        df = CSV.File(files[i]) |> DataFrame
        violin!(ax, fill(i, nrow(df)), mod.(df[:,:μ]))
    end
    #xticks!(ax, xticklabels=["dataset $(i-1)" for i in 1:4])
    return fig, ax
end


function violin_comp(dir1, dir2, label1, label2, quant; mod=(x->1*x))
    dirs1 = filter(isdir, readdir(dir1, join=true))
    files1 = joinpath.(dirs1, basename.(dirs1).*"_$quant.csv")
    dirs2 = filter(isdir, readdir(dir2, join=true))
    files2 = joinpath.(dirs2, basename.(dirs2).*"_$quant.csv")
    fig = Figure(resolution=(1000, 600))
    ax = Axis(fig[1,1],
              xticks = (collect(1:length(files1)), ["dataset $(i-1)" for i in 1:length(files1)]),
              xticklabelrotation=π/2)
    for i in eachindex(files1,files2)
        df1 = CSV.File(files1[i]) |> DataFrame
        df2 = CSV.File(files2[i]) |> DataFrame
        violin!(ax, fill(i, nrow(df1)), mod.(df1[:,:μ]),
                side=:left, color=:blue,
                label=label1)
        violin!(ax, fill(i, nrow(df2)), mod.(df2[:,:μ]),
                side=:right, color=:green,
                label=label2)
    end
    axislegend(ax, [ax.scene.plots[1], ax.scene.plots[2]], [label1, label2])
    #xticks!(ax, xticklabels=["dataset $(i-1)" for i in 1:4])
    return fig, ax
end

function violinmain()
    fig, ax = violin_comp("../RunsCal/SMRing3wf/", "../RunsCal/SMRing2wf/",
                          "M=3", "M=2",
                          "diam")
    ax.ylabel = "diameter (μas)"
    save("../RunsCal/mean_diameter.png", fig)

    fig, ax = violin_comp("../RunsCal/SMRing3wf/", "../RunsCal/SMRing2wf/",
                          "M=3", "M=2",
                          "fwhm")
    ax.ylabel = "width (μas)"
    save("../RunsCal/mean_width.png", fig)

    fig, ax = violin_comp("../RunsCal/SMRing3wf/", "../RunsCal/SMRing2wf/",
                          "M=3", "M=2",
                          "mp1"; mod=(x->180.0/π*x))
    ax.ylabel = "Position Angle (deg)"
    save("../RunsCal/mean_mp1.png", fig)

    fig, ax = violin_comp("../RunsCal/SMRing3wf/", "../RunsCal/SMRing2wf/",
                          "M=3", "M=2",
                          "floor")
    ax.ylabel = "disk flux fraction"
    save("../RunsCal/mean_floor.png", fig)

end
