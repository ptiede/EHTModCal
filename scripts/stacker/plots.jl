using CairoMakie
using ROSESoss
using Measurements

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

function plot_stacker_summary(chain, echain, label, min, max)
    fig = Figure(resolution=(600, 800))
    ax = Axis(fig, ylabel="Probability Density")
    axsa = Axis(fig)
    kμ = kde(echain.μ)
    lower = fill(0.0, length(kμ.x))
    upper = kμ.density
    lines!(ax, kμ.x, upper, color=(:black, 0.2))
    band!(ax, kμ.x, lower, upper, color=(:blue, 0.5), label="μ $label")
    ridge2(axsa, chain)
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

function plotchi2(dir)
    df =  CSV.File(joinpath(dir, "merged_stats.csv")) |> DataFrame

    fig = Figure(reoslution=(300,400))
    ax1 = Axis(fig[1,1], xlabel="Time (UTC)", ylabel="χ² amp.")
    ax2 = Axis(fig[2,1], xlabel="Time (UTC)", ylabel="χ² cp")
    scatter!(ax1, df[!,:time], df[!,:mampchi2])
    scatter!(ax2, df[!,:time], df[!,:mcpchi2])

    linkxaxes!(ax1,ax2)
    hidexdecorations!(ax1, grid=false, ticks=false)

    return fig
end

function create_summary(dir, tmin=0.0, tmax=24.0)
    mdirs = filter(isdir, readdir(dir, join=true))
    println(basename.(mdirs))

    df = DataFrame[]
    for m in mdirs
        mins, maxs, wrapped, quants, labels = parsechainpath(m)
        chain = ChainH5(joinpath(m, "chain.h5"), quants)

    end
end

function StatsBase.mean_and_std(chain::ChainH5)

end


function plotlogz(dir, tavg; tmin=0.0, tmax=1e30, prepend="")
    dirs = filter(isdir, readdir(dir, join=true))
    cfiles =  joinpath.(dirs, Ref("chain.h5"))
    logz = Float64[]
    sims = []
    for (i,c) in pairs(cfiles)
        chain = restricttime(ChainH5(c, (:img_diam, :img_fwhm), 1), tmin, tmax)
        println(length(chain.logz))
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
    #ylims!(ax, -50, 10.0)
    save(joinpath(dir, "logz_"*prepend*".png"), fig)
    return fig
end


function ridge2(ax, chain; amp=3.0, mod=(x->x), kdekwargs=(), )
    times = chain.times
    my = Base.map(mod, chain.chain)
    xlow = minimum(minimum.(my))
    xhigh = maximum(maximum.(my))
    scl = amp*(xhigh-xlow)/360
    kdes = []
    for i in eachindex(times)
        k = kde(my[i]; kdekwargs...)
        push!(kdes, k)
        if k.x[1] < xlow
            xlow = k.x[1]
        end

        if k.x[end] > xhigh
            xhigh = k.x[end]
        end
    end
    band!(ax, [xlow,xhigh], fill(13.25,2), fill(13.7,2), color=(:grey, 0.4))
    for t in times
        lines!(ax, [xlow, xhigh], fill(t,2), color=(:grey, 0.5), order=-1)
    end


    for i in reverse(eachindex(times))
        k = kdes[i]
        lower = fill(times[i], length(k.x))
        upper = times[i] .+ k.density*scl
        imax = argmax(upper)
        band!(ax, k.x, lower, upper, linewidth=2.0, transparency=true, color=k.density, colormap=(:viridis, 0.75))
        scatter!(ax, [k.x[imax]], [upper[imax]], color=:red, markersize=5.0)
        lines!(ax, k.x, upper, linewidth=0.5, color=(:black, 0.8))
    end

    xlims!(xlow, xhigh)
    ax.yticks = times[1:4:end]
    ax.ylabel = "Time (UT)"
    return ax
end

function ridge2(chain::ChainH5; amp=3.0, mod=(x->x), kdekwargs=(),)
    fig = Figure(resolution=(600,600))
    ax = Axis(fig[1,1])
    ax = ridge2(ax, chain; amp, mod, kdekwargs)
    return fig
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

function ridge(chain::ChainH5; mod=(x->x), kdekwargs=(),)
    fig = Figure(resolution=(700,700))
    ax = Axis(fig[1,1])
    ax = ridge(ax, chain, ;mod, kdekwargs)
    return fig
end



function plotmeanim(model, tv, nsamples)
    mim = meanframe(model, tv, nsamples)
    fig = Figure(resolution=(400,400))
    ax = Axis(fig[1,1], aspect=DataAspect(),xlabel="RA (μas)", ylabel="DEC (μas)", xreversed=true)
    image!(ax, ROSE.imagepixels(mim)..., mim', colormap=Reverse(:lajolla))
    return fig
end

function plotmeanim(ax, model, tv, nsamples)
    mim = meanframe(model, tv, nsamples)
    image!(ax, ROSE.imagepixels(mim)..., mim', colormap=Reverse(:lajolla))
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

function generatemovie(dir; tmin=12.6, tmax=14.1, nsamples=100, prepend="")
    model = loadmodel(dir)
    println(model)
    quants = parsechainpath(dir)[4]
    chain = ChainH5(joinpath(dir, "chain.h5"), quants, 1)
    chain = restricttime(chain, 12.6, 14.1)
    meanmovie(model, chain, nsamples, joinpath(dir, "mean_movie_$prepend.gif"))
end

function plotcollage(m, tv, dims=(3,3))
    fig = Figure(resolution=(800,800))
    axes = [Axis(fig[i,j], xreversed=true, aspect=DataAspect()) for i in 1:dims[1], j in 1:dims[2]]
    for ax in axes
        ind = rand(1:length(tv))
        plotim!(ax, m, tv[ind])
        hidedecorations!(ax)
    end
    return fig
end

function summarize_sgra(dir; tmin=0.0, tmax=24.0, genmovie=false, prepend="")
    dirs = filter(isdir, readdir(dir, join=true))
    cfiles =  joinpath.(dirs, Ref("chain.h5"))
    for c in cfiles
        summarize_snapshot(c; tmin, tmax, genmovie, prepend)
    end

end



function summarize_snapshot(cfile; tmin=0.0, tmax=24.0, genmovie=false, prepend="")
    mins, maxs, wrapped, quants, labels = parsechainpath(dirname(cfile))
    chain = restricttime(ChainH5(cfile, quants), tmin, tmax)
    #println(chain)
    for i in eachindex(quants)
        println(quants[i])
        if occursin("img_mp_", String(quants[i]))
            fig = ridge2(getparam(chain, quants[i]), mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180,180),))
            fig.content[1].xticks=-150:50:150
            xlims!(-180, 180)
        else
            fig = ridge2(getparam(chain, quants[i]))
        end
        fig.content[1].xlabel = labels[i]
        display(fig)
        save(joinpath(dirname(cfile), "ridge_$(String(quants[i]))_$prepend.png"), fig)
    end
    if genmovie
        generatemovie(dirname(cfile); tmin, tmax, prepend)
    end
end

function summarize_ha(cfile)
    chainha = CSV.File(cfile) |> DataFrame
    tv = df2tv(chainha)
    m = loadmodel(dirname(dirname(cfile)))
    println(m)

    fig = plotmeanim(m, tv, 100)
    fig.current_axis.x.title = basename(dirname(dirname(cfile)))
    display(fig)
    save(joinpath(dirname(cfile), "mean_image.png"), fig)

    fig = plotcollage(m, tv)
    display(fig)
    save(joinpath(dirname(cfile), "collage.png"), fig)


    fig = Figure(resolution=(400,400))
    ax = Axis(fig[1,1], xlabel="Diameter", ylabel="fractional width")
    scatter!(ax, tv.diam, tv.fwhm./tv.diam, label="Delta diam")
    diamdb = @. tv.diam - 1/(4*log(2))*tv.fwhm^2/tv.diam
    scatter!(ax, diamdb, tv.fwhm./tv.diam, label="Peak diam")

    if occursin("_floor_order", cfile)
        scatter!(ax, tv.diam./(1 .+ 0.6*tv.floor), tv.fwhm./tv.diam, label="Effective diam")
    end
    axislegend(ax)
    xlims!(ax, 30.0, 70.0)
    ylims!(ax, 0.0, 0.8)
    display(fig)
    save(joinpath(dirname(cfile), "corr_image.png"), fig)

    #fig = plotchi2(dirname(dirname(cfile)))
    #display(fig)
    #save(joinpath(dirname(cfile), "chi2_image.png"), fig)

end

function plotim(model, p::NamedTuple)
    img = Soss.predict(model, p)
    mim = intensitymap(img, 256, 256, 150.0, 150.0)
    fig = Figure(resolution=(400,400))
    ax = Axis(fig[1,1], aspect=DataAspect(),xlabel="RA (μas)", ylabel="DEC (μas)", xreversed=true)
    image!(ax, ROSE.imagepixels(mim)..., mim', colormap=Reverse(:lajolla))
    return fig
end


function plotim!(ax, model, p::NamedTuple)
    img = Soss.predict(model, p)
    mim = intensitymap(img, 256, 256, 150.0, 150.0)
    image!(ax, ROSE.imagepixels(mim)..., mim', colormap=Reverse(:lajolla))
    lines!(ax, [10,60],[-60, -60], color=:white)
    text!(ax, "50 μas", position = (35, -55), align=(:center, :baseline), color=:white, textsize=14)
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
