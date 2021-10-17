using CairoMakie
using EHTModelStacker
using KernelDensity

function ridge2(ax, chain; amp=3.0, mod=(x->x), kdekwargs=(), color=(RGB(0.0,0.4,0.6), 0.5))
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
    band!(ax, [xlow,xhigh], fill(13.25,2), fill(13.7,2), color=(:grey, 0.5))
    #for t in times
    #    lines!(ax, [xlow, xhigh], fill(t,2), color=(:grey, 0.3), order=-1)
    #end


    for i in reverse(eachindex(times))
        k = kdes[i]
        lower = fill(times[i], length(k.x))
        upper = times[i] .+ k.density*scl
        imax = argmax(upper)
        band!(ax, k.x, lower, upper, linewidth=0.1, color=color)
        lines!(ax, k.x, upper, linewidth=0.5, color=(:black, 0.6))
    end

    xlims!(xlow, xhigh)
    return ax
end

function generateticks(locs)
    ticklbls = [L"$ %$x $" for x in locs]
    return locs, ticklbls
end


models = filter(isdir, readdir("_research/SgrARunsFinal_BT/hops/3599/LO/snapshot_60/noisefrac0.02/", join=true))
c1_99 = restricttime(ChainH5(joinpath(models[1], "chain.h5"), (:img_diam, :img_mp_1)), 12.455, 14.2)
c2_99 = restricttime(ChainH5(joinpath(models[2], "chain.h5"), (:img_diam, :img_mp_1)), 12.455, 14.2)
c3_99 = restricttime(ChainH5(joinpath(models[3], "chain.h5"), (:img_diam, :img_mp_1)), 12.455, 14.2)
c4_99 = restricttime(ChainH5(joinpath(models[4], "chain.h5"), (:img_diam, :img_mp_1)), 12.455, 14.2)

models = filter(isdir, readdir("_research/SgrARunsFinal_BT/hops/3598/LO/snapshot_60/noisefrac0.02/", join=true))
c1_98 = restricttime(ChainH5(joinpath(models[1], "chain.h5"), (:img_diam, :img_mp_1)), 12.455, 14.2)
c2_98 = restricttime(ChainH5(joinpath(models[2], "chain.h5"), (:img_diam, :img_mp_1)), 12.455, 14.2)
c3_98 = restricttime(ChainH5(joinpath(models[3], "chain.h5"), (:img_diam, :img_mp_1)), 12.455, 14.2)
c4_98 = restricttime(ChainH5(joinpath(models[4], "chain.h5"), (:img_diam, :img_mp_1)), 12.455, 14.2)

logzs_98 = sum.([c1_98.logz, c2_98.logz, c3_98.logz, c4_98.logz])
dlogz_98 = -maximum(logzs_98) .+ logzs_98
logzs_99 = sum.([c1_99.logz, c2_99.logz, c3_99.logz, c4_99.logz])
dlogz_99 = -maximum(logzs_99) .+ logzs_99


color98 = RGB(0.0,0.4,0.6)
color99 = RGB(0.8, 0.5, 0.0)

fig = Figure(resolution=(1000,800))
h = fig[1,1] = GridLayout()
g = fig[2,1] = GridLayout()

rowsize!(fig.layout, 1, Relative(1/5))
yticklbls = [L"$ %$x $" for x in -300:100:25]
axt = Axis(h[1,1], xticks=(1:4, [L"m=1", L"m=2", L"m=3", L"m=4"]),
                   yticks=generateticks(-300:100:25),
                   ylabel = L"$\Delta$log$Z$",
                   limits=((0.5,4.5), (-300,25)))


scatter!(axt, [1,2,3,4], dlogz_98; color=color98)
scatter!(axt, [1,2,3,4], dlogz_99; color=color99)


#hidedecorations!(axt1, label=true, ticklabels=true, ticks=true, grid=true)
xtickloc = -150:50:150
xticklbls = [L"$ %$x^\circ$" for x in xtickloc]
ax1 = Axis(g[1,1], xticks=xtickloc,
                   yticks = generateticks(c1.times[1:6:end]), limits=((-180,180), (c1.times[1], 14.4)),
                   ylabel = L"April $6$ (UT)", ygridvisible=false
                   )
ridge2(ax1, getparam(c1_98, :img_mp_1), amp=2.5, mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180.0,180.0),))
hidexdecorations!(ax1, ticklabels=true, ticks=false, grid=false)


ax2 = Axis(g[1,2], xticks=xtickloc,
                   limits=((-180,180), (c1.times[1], 14.4)))
ridge2(ax2, getparam(c2_98, :img_mp_1),amp=2.5, mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180.0,180.0),))
hidexdecorations!(ax2, ticklabels=true, ticks=false, grid=false)
hideydecorations!(ax2, ticklabels=true, ticks=false)

ax3 = Axis(g[1,3], xticks=xtickloc,
                    limits=((-180,180), (c1.times[1], 14.4)))
ridge2(ax3, getparam(c3_98, :img_mp_1), mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180.0,180.0),))
hidexdecorations!(ax3, ticklabels=true, ticks=false, grid=false)
hideydecorations!(ax3, ticklabels=true, ticks=false)

ax4 = Axis(g[1,4], xticks=xtickloc,
                    limits=((-180,180), (c1.times[1], 14.4)))
ridge2(ax4, getparam(c4_98, :img_mp_1), amp=4.0, mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180.0,180.0),))
hidexdecorations!(ax4, ticklabels=true, ticks=false, grid=false)
hideydecorations!(ax4, ticklabels=true, ticks=false)

ax1 = Axis(g[2,1], xticks=(xtickloc, xticklbls),
                   yticks = generateticks(c1.times[1:6:end]), limits=((-180,180), (c1.times[1], 14.4)),
                   ylabel=L"April $7$ (UT)",  xticklabelrotation=π/4,  ygridvisibile=false)
ridge2(ax1, getparam(c1_99, :img_mp_1), mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180.0,180.0),); color=(color99, :0.5))

ax2 = Axis(g[2,2], xticks=(xtickloc, xticklbls), xticklabelrotation=π/4, ygridvisibile=false,
                   limits=((-180,180), (c1.times[1], 14.4)))
ridge2(ax2, getparam(c2_99, :img_mp_1), amp=4.0, mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180.0,180.0),); color=(color99, :0.5))
hideydecorations!(ax2, label=true, ticklabels=true, ticks=false)

ax3 = Axis(g[2,3], xticks=(xtickloc, xticklbls), xticklabelrotation=π/4, ygridvisibile=false,
                    limits=((-180,180), (c1.times[1], 14.4)))
ridge2(ax3, getparam(c3_99, :img_mp_1), amp=5.0, mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180.0,180.0),);color=(color99, :0.5))
hideydecorations!(ax3, label=true, ticklabels=true, ticks=false)

ax4 = Axis(g[2,4], xticks=(xtickloc, xticklbls), xticklabelrotation=π/4,  ygridvisibile=false,
                    limits=((-180,180), (c1.times[1], 14.4)))
ridge2(ax4, getparam(c4_99, :img_mp_1), amp=5.0, mod=(x->-rad2deg.(x)), kdekwargs=(boundary=(-180.0,180.0),); color=(color99, :0.5))
hideydecorations!(ax4, label=true, ticklabels=true, ticks=false)
#linkaxes!(ax1,ax2,ax3,ax4)


rowgap!(fig.layout, 4.0)
Label(g[end+1,:], L"Position Angle $E$ of $N$")
rowgap!(g, 1.0)
colgap!(g, 10.0)
fig
save("model_summary.pdf", fig)
