using CairoMakie

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

function plotlogz(dir, tavg; tmin=0.0, tmax=1e30)
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
    save(joinpath(dir, "logz.png"), fig)
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

function summarize_ha(cfile)
    chainha = CSV.File(cfile) |> DataFrame
    tv = df2tv(chainha)
    m = loadmodel(dirname(dirname(cfile)))

    fig = plotmeanim(m, tv, 100)
    fig.current_axis.x.title = basename(dirname(dirname(cfile)))
    display(fig)
    save(joinpath(dirname(cfile), "mean_image.png"), fig)

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
    image!(ax, ROSE.imagepixels(mim)..., mim', colormap=:afmhot)
    return fig
end




function comboavgplot(fchainsa, fcsa, outname)
    chainsa = ChainH5(fchainsa, "diam", 2000)
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
