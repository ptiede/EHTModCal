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



function build_model_all(θ::NamedTuple, mins, maxs)
    μ, Σ = θ.μ, θ.σ
    transition = EHTModelStacker.MvNormalFast(μ, Σ.^2)
    prior = EHTModelStacker.MvUniform(mins, maxs)
    return SnapshotWeights(transition, prior)
end


function build_loglklhd(build_model, chain, min, max)
    return function (θ::NamedTuple)
        ws = build_model(θ, min, max)
        return lpdf(ws, chain)
    end
end

function generatemovie(dir; bt=true, nsamples=100)
    model = loadmodel(dir)
    println(model)
    quants = parsechainpath(dir)[3]
    chain = ChainH5(joinpath(dir, "chain.h5"), quants, 1)
    if bt
        chain = restricttime(chain, 12.6, 14.1)
    end
    meanmovie(model, chain, nsamples, joinpath(dir, "mean_movie.gif"))
end

function getsnapshot_tv(chain::ChainH5, i::Int)
    k = keys(chain)
    nms = String.(k)
    iimg = findall(x->occursin("img_",x), nms)
    inon = findall(x->!((occursin("img_mp_",x)||occursin("img_ma_",x))&&occursin("img_", x)), nms)
    order = length(findall(x->occursin("img_ma_",x), nms))
    ks = Symbol.(last.(split.(nms[inon], Ref("img_"))))
    ntproto = namedtuple(ks..., :ma, :mp)
    ma = nestedview(getparam(chain, Tuple(Symbol.(["img_ma_$i" for i in 1:order]))).chain[i])
    mp = nestedview(getparam(chain, Tuple(Symbol.(["img_mp_$i" for i in 1:order]))).chain[i])
    imgt = Tuple(chain.chain[i][k,:] for k in inon)
    return TupleVector(ntproto((imgt..., ma, mp)))
end


function loadmodel(dir)
    model = split(dir, "/")[end]
    order = parse(Int, split(model, "-")[end])

    if occursin("mring_floor_order", model)
        img = ROSESoss.mringwfloor(N=order,)
    elseif occursin("mring_order", model)
        img = ROSESoss.mring(N=order,)
    elseif occursin("mring_floor_stretch_order", model)
        img = ROSESoss.smringwfloor(N=order,)
    else
        throw("model $model not found")
    end
    return img
end

function plotmeanim(model, tv, nsamples)
    mim = meanframe(model, tv, nsamples)
    fig = Figure(resolution=(400,400))
    ax = Axis(fig[1,1], aspect=DataAspect(),xlabel="RA (μas)", ylabel="DEC (μas)", xreversed=true)
    image!(ax, ROSE.imagepixels(mim)..., mim', colormap=:afmhot)
    return fig
end

function plotmeanim(ax, model, tv, nsamples)
    mim = meanframe(model, tv, nsamples)
    image!(ax, ROSE.imagepixels(mim)..., mim', colormap=:afmhot)
end


function meanmovie(model, chain, nsamples, outname; framerate=5)
    fig = Figure(resolution=(400,400))
    ax = Axis(fig[1,1], xlabel="RA (μas)", ylabel="DEC (μas)", xreversed=true)
    record(fig, outname, eachindex(chain.times); framerate) do i
        tv = getsnapshot_tv(chain, i)
        ax.title="Time = $(chain.times[i])"
        plotmeanim(ax, model, tv, nsamples)
    end
end

function meanframe(model, chain::ChainH5, sindx::Int, nsamples::Int)
    tv = getsnapshot_tv(chain, sindx)
    return meanframe(model, tv, nsamples)
end

function meanframe(model, tv::TupleVector, nsamples::Int; fov=150, npix=256)
    inds = rand(1:length(tv), nsamples)
    simacc = ROSE.StokesImage(zeros(npix, npix), fov, fov)
    tmp = similar(simacc)
    tst = tv[1]
    for i in inds
        make_sim!(tmp, model, tst)
        simacc += tmp
    end
    return simacc./nsamples
end

function make_sim!(tmp, model, p::NamedTuple; fov=150, npix=256)
    img = Soss.predict(model, p)
    intensitymap!(tmp, img)
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


function average_chain_diag(build_model, mins, maxs, chain; nlive=1500)
    priors = (μ = Product(Uniform.(mins, maxs)), σ = Product(Uniform.(0.0, maxs .- mins)))
    t = ascube(priors)
    p0 = HypercubeTransform.transform(t, fill(0.5, dimension(t)))
    fp0, unflatten = ParameterHandling.flatten(p0)
    lp = build_loglklhd(build_model, chain, mins, maxs)∘unflatten
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



function plot_stacker_summary(chain, echain, label, min, max)
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
    xlims!(ax, min, max)
    return fig
end

function _mkdf(echain, keys)
    df = DataFrame()
    for (i,k) in enumerate(keys)
        insertcols!(df, Symbol("μ_"*String(k))=>getindex.(echain.μ, i))
    end
    for (i,k) in enumerate(keys)
        insertcols!(df, Symbol("σ_"*String(k))=>getindex.(echain.σ, i))
    end
    df
end


function parsechainpath(cfile)
    model = split(cfile, "/")[end]
    mins = [25.0, #diam
            1.0, #width
            0.0, #flux
            ]
    maxs = [85.0,
            40.0,
            5.0
           ]
    quants = [:img_diam, :img_fwhm, :img_f]
    labels = ["diameter (μas)", "width (μas)", "flux (Jy)"]

    if occursin("gfloor", model)
        push!(mins, 0.0, 10.0)
        push!(maxs, 1.0, 90.0)
        push!(quants, :img_floor, :img_dg)
        push!(labels, "floor flux fraction", "Gaussian diameter")
    elseif occursin("gfloor", model)
        push!(mins, 0.0)
        push!(maxs, 1.0)
        push!(quants, :img_floor)
        push!(labels, "floor flux fraction")
    end

    if occursin("stretch", model)
        push!(mins, 0.0, π/2)
        push!(maxs, 0.5, -π/2)
        push!(quants, :img_τ, :img_ξτ)
        push!(labels, "ellipticity τ", "ellipticity PA (rad) W of N")
    end

    order = parse(Int, split(model, "-")[end])

    for i in 1:order
        push!(quants, Symbol("img_ma_$i"), Symbol("img_mp_$i"))
        push!(mins, 0.0, -1π)
        push!(maxs, 0.5, 1π)
        push!(labels, "m=$i amp", "m=1 phase (rad) W of N")
    end

    return mins, maxs, Tuple(quants), labels

end




function quickprocess(cfile; plotres=true, tmin=0, tmax=1e30)
    mins, maxs, quants, labels = parsechainpath(cfile)
    outdir = joinpath(cfile, "ChainHA")

    process(joinpath(cfile, "chain.h5"), mins, maxs, quants, labels, outdir; plotres, tmin, tmax)

end


function process(cfile, mins, maxs, quants, labels, outdir; plotres=true, tmin=0, tmax=1e30)
    #diameter
    mkpath(outdir)
    chain = ChainH5(cfile, quants, 500)
    chainall = restricttime(chain, tmin, tmax)
    echain = average_chain_diag(build_model_all, mins, maxs, chainall)
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

function ridge(ax, chain::ChainH5; mod=(x->x), kdekwargs=(),)
    times = chain.times
    height = (times[end]-times[1])/length(times)*2.0
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

function plotlogz(dir, tavg; tmin=0.0, tmax=1e30)
    dirs = filter(isdir, readdir(dir, join=true))
    cfiles =  joinpath.(dirs, Ref("chain.h5"))
    logz = Float64[]
    sims = []
    for (i,c) in pairs(cfiles)
        chain = restricttime(ChainH5(c, (:img_diam, :img_fwhm), 1), tmin, tmax)
        push!(logz, sum(chain.logz))
        push!(sims, basename(dirs[i]))
    end
    fig = Figure(resolution=(800,800))
    ax = Axis(fig[1,1],
            xticks = (collect(1:length(logz)), sims),
            xticklabelrotation=π/4,
            ylabel="Δlogz",
            title=tavg)
    imax = argmax(logz)
    scatter!(ax,1:length(logz), logz .- maximum(logz))
    scatter!(ax, imax:imax, [0.0], color=:orange)
    ylims!(ax, -50, 10.0)
    save(joinpath(dir, "logz.png"), fig)
    return fig

end


function comboavgplot(fchainsa, fcsa, outname)
    chainsa = ChainH5(fchainsa, "img_diam")
    csa = CSV.File(fcsa) |> DataFrame

    fig = Figure(resolution=(1200, 1600))
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
