using ROSE: StokesImage
using CSV, DataFrames
using Statistics
using KernelDensity
using ROSE
using Soss
using CairoMakie
using Printf
using PyCall
using StatsBase
using HDF5


const fwhmfac = 2*sqrt(2*log(2))
function parsemodel(pfile)
    data, model = open(pfile, "r") do io
        data = readline(io)
        readline(io)
        model = readuntil(io, "\nfit_gains")
        return data, model
    end
    dfile = last(split(data, ": "))
    return dfile,  eval(Meta.parse(model))
end




function hdf5_chain(dir, outname)
    dfchain, dfsum = load_chains(dir)
    h5open(outname, "w") do fid
        fid["time"] = dfsum[:,:time]
        fid["logz"] = dfsum[:,:logz]
        pid = create_group(fid, "params")
        for i in 1:length(dfchain)
            sid = create_group(pid, "scan$i")
            keys = names(dfchain[i])
            for k in keys
                write(sid, k, dfchain[i][:,Symbol(k)])
            end
        end
    end
end


function extract_params(file, outdir, quant; nsamples=2000)
    open(joinpath(outdir, "$(quant).dat"), "w") do io
        h5open(file, "r") do fid
            times = read(fid["time"])
            params = fid["params"]
            for i in 1:length(times)
                p = read(params["scan$i"][quant])[1:nsamples]
                writedlm(io, p', " ")
            end
        end
    end
end

function find_bhm(simdir)
    sfile = "diamdb.dat_out_gauss/info/post_summary.csv"
    dirs = filter(isdir, readdir(simdir, join=true))
    df = DataFrame()
    for d in dirs
        sdf = CSV.File(joinpath(d, sfile)) |> DataFrame
        append!(df, sdf)
    end
    return df
end

function mainextract(simdir)
    hfiles = filter(endswith(".h5"),readdir(simdir, join=true))
    for h in hfiles
        od = first(split(h, "_"))
        extract_params(h, od, "diamdb")
        hdir = pwd()
        try
            cd(od)
            run(`posteriorstacker.py diamdb.dat 10 80 10 --name diameter`)
            run(`convert -density 100 diamdb.dat_out.pdf diamdb.png`)
        finally
            cd(hdir)
        end
    end
end


function plotframe!(ax, mim, stats)

    image!(ax, imagepixels(mim)..., mim', colormap=:viridis)
    ax.title="Time = $(stats[1,:time])"
    #text!(ax, srampchi2, position=(70.0, 65.0), color=:white, textsize=12)
    #text!(ax, smampchi2, position=(70.0, 50.0), color=:white, textsize=12)
    #text!(ax, srcpchi2, position=(70.0, -55.0), color=:white, textsize=12)
    #text!(ax, smcpchi2, position=(70.0, -70.0), color=:white, textsize=12)
    #text!(ax, snbl, position=(50, 60), color=:white, textsize=10)
    hidedecorations!(ax)
    return ax
end

function collage(dir; nimages = 5)
    aratio = 1.0
    fig = Figure(resolution=(aratio*500*2.5,500*2.5))
    cfiles = filter(endswith(".csv"),readdir(joinpath(dir,"Chain"), join=true))
    sfiles = filter(endswith(".csv"),readdir(joinpath(dir, "Stats"), join=true))
    #Now because I indexed these in a dumb way lets reorder everything
    ind =  parse.(Int, last.(split.(first.(splitext.(basename.(cfiles))), "scan-")))
    sind = sortperm(ind)
    ny = Int(ceil(sqrt(1/aratio*length(sind))))
    nx = Int(ceil(aratio*ny))
    cifs = CartesianIndices((1:nx,1:ny))
    dfile, m = parsemodel(joinpath(dir, "params-0.dat"))
    pr = prior(m, :amp, :cphase, :g)
    primg(x) = getfield(Soss.predict(pr, x), :img)
    for (i,is) in enumerate(sind)
        dd = replace(dfile, "scan0"=>"scan$(i-1)")
        obs = ehtim.obsdata.load_uvfits(dd)
        nbl = length(obs.bllist())
        mim = ROSE.StokesImage(zeros(128,128), 160.0, 160.0)
        fim = deepcopy(mim)
        dfchain = CSV.File(cfiles[is]) |> DataFrame
        dfstats = CSV.File(sfiles[is]) |> DataFrame
        #dfstats.nbl = [-1]
        for n in 1:nimages
            model = primg(copy(dfchain[n,3:(end-1)]))
            mim += intensitymap!(fim, model)
        end
        mim = mim./nimages
        ix,iy = Tuple(cifs[i])
        ax = Axis(fig[iy,ix], xreversed=true, aspect=DataAspect())
        plotframe!(ax, mim, dfstats)
    end
    CairoMakie.trim!(fig.layout)
    return fig
end



function collage(bims::AbstractVector, dfstats)
    aratio = Base.MathConstants.golden
    fig = Figure(resolution=(aratio*500*2.5,500*2.5))
    ny = Int(ceil(sqrt(1/aratio*length(bims))))
    nx = Int(ceil(aratio*ny))
    cifs = CartesianIndices((1:nx,1:ny))
    for (i,mim) in enumerate(bims)
        ix,iy = Tuple(cifs[i])
        dfstats.nbl = zeros(nrow(dfstats))
        plotframe!(fig[iy,ix], mim, dfstats)
    end
    return fig
end




function merge_summaries(sumdir)
    files = filter(endswith(".csv"), readdir(sumdir, join=true))
    df = CSV.File(files[1]) |> DataFrame
    for f in files[2:end]
        append!(df, CSV.File(f)|> DataFrame)
    end
    return sort!(df, :time)
end


function plot_res(dir; fitgains=true, pcollage=true)
    dfchains, dfsum = load_chains(dir)
    if pcollage
        fig = collage(dir; nimages=250)
        save(joinpath(dir, "collage.png"), fig)
    end

    if fitgains
       fig = plot_mring_gains(dfsum, dfchains)
       save(joinpath(dir, "gains.png"), fig)
    end

    fig = Figure(resolution=(500,500/1.5))
    ax = Axis(fig[1,1], ylabel = "log Z", xlabel = "Time UTC")
    CairoMakie.errorbars!(ax, dfsum[:,:time], dfsum[:,:logz], dfsum[:,:logzerr], whiskerwidth=1, linewidth=3)
    CairoMakie.scatter!(ax, dfsum[:,:time], dfsum[:,:logz])
    save(joinpath(dir, "logz.png"), fig)

    fig = Figure(resolution=(500,500/1.5))
    ax = Axis(fig[1,1], ylabel = "reduced χ²", xlabel = "Time UTC")
    CairoMakie.scatter!(ax, dfsum[:,:time], dfsum[:,:rchi2], color=:green, label="Total")
    CairoMakie.scatter!(ax, dfsum[:,:time], dfsum[:,:rampchi2], color=:blue, label="Amp")
    CairoMakie.scatter!(ax, dfsum[:,:time], dfsum[:,:rcpchi2], color=:orange, label="CP")
    CairoMakie.ylims!(ax,(0.0,10.0))
    axislegend(ax, framevisible=false, position=:lt)
    save(joinpath(dir, "rchi2.png"), fig)

    fig = Figure(resolution=(500,500/1.5))
    ax = Axis(fig[1,1], ylabel = "mean χ²", xlabel = "Time UTC")
    CairoMakie.scatter!(ax, dfsum[:,:time], dfsum[:,:mampchi2], color=:blue, label="Amp")
    CairoMakie.scatter!(ax, dfsum[:,:time], dfsum[:,:mcpchi2], color=:orange, label="CP")
    axislegend(ax, framevisible=false, position=:lt)
    save(joinpath(dir, "mchi2.png"), fig)


    fig = Figure(resolution=(750,800))
    _,ax = ridge(fig[1,1],dfsum[:,:time], dfchains, :diamdb)
    CairoMakie.xlims!(ax, 38.0, 72.0)
    ax.xlabel = "blurred Diameter [μas]"

    _,ax = ridge(fig[1,2],dfsum[:,:time], dfchains, :diam)
    CairoMakie.xlims!(ax, 38.0, 72.0)
    ax.xlabel = "delta Diameter [μas]"

    _,ax = ridge(fig[1,3], dfsum[:,:time], dfchains, :fwhm)
    ax.xlabel = "FWHM [μas]"

    _,ax = ridge(fig[2,1], dfsum[:,:time], dfchains, :mp1, mod=(x->-x*180.0/π), kdekwargs=(boundary=(-180.0,180.0),))
    ax.xlabel = "Position Angle [deg]"

    _,ax = ridge(fig[2,2], dfsum[:,:time], dfchains, :ma1, mod=(x->x*2))
    CairoMakie.xlims!(ax, (0.0, 1.0))
    ax.xlabel = "Slash"

    try
        _,ax = ridge(fig[2,3], dfsum[:,:time], dfchains, :floor)
        CairoMakie.xlims!(ax, (0.0, 1.0))
        ax.xlabel = "Floor"
    catch
        @info "Model has no floor, skipping"
    end

    try
        _,ax = ridge(fig[3,1], dfsum[:,:time], dfchains, :τ)
        CairoMakie.xlims!(ax, (0.0, 1.0))
        ax.xlabel = "ellipticity τ"
        _,ax = ridge(fig[3,2], dfsum[:,:time], dfchains, :ξτ, mod=(x->-x*180.0/π), kdekwargs=(boundary=(-90.0,90.0),))
        #CairoMakie.xlims!(ax, (0.0, 1.0))
        ax.xlabel = "ellipticity angle ξτ"
        #display(fig)
    catch
        @info "Model has not ellipticity skipping"
    end



    #ax = Axis(fig[3,3], ylabel = "reduced χ²", xlabel = "Time UTC")
    #CairoMakie.scatter!(ax, dfsum[:,:time], dfsum[:,:rampchi2], color=:blue, label="Amp")
    #CairoMakie.scatter!(ax, dfsum[:,:time], dfsum[:,:rcpchi2], color=:orange, label="CP")
    #axislegend(ax, framevisible=false, position=:lt)
    save(joinpath(dir, "params.png"), fig)

end


function ridge_quant(fig, times, qchain, xlabel; mod=(x->x), kdekwargs=(),)
    height = (times[end]-times[1])/length(times)*1.5
    y = qchain
    my = map(mod, y)
    ax = Axis(fig)
    cls = to_colormap(:viridis, length(times))
    for i in 1:length(times)
        k = kde(my[i]; kdekwargs...)
        lower = fill(times[i], length(k.x))
        upper = times[i] .+ k.density*height/maximum(k.density)
        lines!(k.x, upper, color=(:black, 0.9))
        band!(k.x, lower, upper, color=(cls[i], 0.5))
    end
    ax.ylabel = "Time UTC"
    ax.xlabel = xlabel
    #ax.yticks = (range(times[1],times[end-1], length=10))
    return fig, ax
end


function ridge(fig, times, dfchain, quant; mod=(x->x), kdekwargs=(),)
    height = (times[end]-times[1])/length(times)*3.0
    y = getproperty.(dfchain, quant)
    my = map(mod, y)
    ax = Axis(fig)
    cls = to_colormap(:viridis, length(times))
    for i in 1:length(times)
        k = kde(my[i]; kdekwargs...)
        lower = fill(times[i], length(k.x))
        upper = times[i] .+ k.density*height/maximum(k.density)
        lines!(k.x, upper, color=(:black, 0.2))
        band!(k.x, lower, upper, color=(cls[i], 0.5))
    end
    ax.ylabel = "Time UTC"
    ax.xlabel = "$(String(quant))"
    #ax.yticks = (range(times[1],times[end-1], length=10))
    return fig, ax
end

function ridge!(ax, times, dfchain, quant; mod=(x->x), kdekwargs=(),)
    height = (times[end]-times[1])/length(times)*1.5
    y = getproperty.(dfchain, quant)
    my = map(mod, y)
    cls = to_colormap(:viridis, length(times))
    for i in 1:length(times)
        k = kde(my[i]; kdekwargs...)
        lower = fill(times[i], length(k.x))
        upper = times[i] .+ k.density*height/maximum(k.density)
        #lines!(k.x, upper, color=(:black, 0.9))
        band!(k.x, lower, upper, color=(cls[i], 0.5))
    end
    return ax
end



"""
    load_chains(dir)
Loads the set of csv chains and summary file assuming the structure present in parallel_main
"""
function load_chains(dir)
    files = filter(endswith(".csv"), readdir(joinpath(dir,"Chain"), join=true))
    sfile = joinpath(dir,"merged_stats.csv")
    dfsum = merge_summaries(joinpath(dir, "Stats"))
    CSV.write(sfile, dfsum)

    ind = parse.(Int,first.(splitext.(last.(split.(files, '-')))))
    sind = sortperm(ind)
    dfs = CSV.File.(files) .|> DataFrame
    for df in dfs
        df.diamdb = @. df.diam - 1/(4*log(2))*df.fwhm^2/df.diam
    end

    return dfs[sind], dfsum
end




function plot_mring_gains(dfsum, dfchain)
    ps = []
    stations = [:AP, :AZ, :JC, :SM, :AA, :LM, :SP, :PV]
    fig = Figure(resolution=(1200, 500))
    cifs = CartesianIndices((1:2, 1:4))
    for (i,s) in enumerate(stations)
        iy,ix = Tuple(cifs[i])
        ax = Axis(fig[iy,ix],
                  xlabel="Time UTC",
                  ylabel="gain amp",
                  title=String(s)
                )
        p = Plots.plot()
        mrad = quantile.(getproperty.(dfchain, s), 0.5)
        qup = quantile.(getproperty.(dfchain, s), 0.86)
        qlo = quantile.(getproperty.(dfchain, s), 0.14)
        rangebars!(dfsum[!,:time], qlo, qup)
        CairoMakie.scatter!(dfsum[!,:time], mrad)
    end
    linkaxes!(contents(fig.layout)...)
    hideydecorations!.(contents(fig.layout)[3:end], grid = false, ticks=false)
    return fig
end


function main(sim="RunsSMRing2wf" ,np = "preprocessed_test/scan-average/noisefrac0.02/")
    dirs = filter(isdir, readdir(sim, join=true))
    for d in dirs
        @info "On $d"
        plot_res(joinpath(d, np))
    end
end


function foo()
    dfchainsm2, dfsumm2 = load_chains("HeDataSgrAMring2VACP/")
    dfchainsm2wf, dfsumm2wf = load_chains("HeDataSgrAMring2wfVACP/")
    dfchainssm2wf, dfsumsm2wf = load_chains("HeDataSgrASMring2wfVACP/")

    fig = Figure(resolution=(600,400))
    ax = Axis(fig[1,1], ylabel = "log Z", xlabel = "Time UTC")
    CairoMakie.scatter!(ax, dfsumm2[:, :time],dfsumm2[:, :logz], color=:blue, label="2-ring")
    CairoMakie.scatter!(ax, dfsumm2wf[:, :time],dfsumm2wf[:, :logz], color=:orange, label="2-ring w floor")
    CairoMakie.scatter!(ax, dfsumsm2wf[:, :time],dfsumsm2wf[:, :logz], color=:purple, label="2-ring w floor & stretch")
    axislegend(ax, position=:rt)

    return fig
end

function bmaall_quant(quant, dirs...)
    dfchains = Vector{DataFrame}[]
    dfsums = DataFrame[]
    for dir in dirs
        dfchain, dfsum = load_chains(dir)
        push!(dfchains, dfchain)
        push!(dfsums, dfsum)
    end
    quant_bma = Vector{Float64}[]
    # i is over times
    # j is over models
    for i in 1:length(dfchains[1])
        println(i)
        sums = [dfsums[j][i,:] for j in 1:length(dfsums)]
        quantchain = [dfchains[j][i][:,quant] for j in 1:length(dfchains)]
        cchain, weights = bma_ss_quant(quantchain, sums, [dirs...])
        println(weights[1])
        echain = sample(cchain, Weights(weights), 3000)
        push!(quant_bma, echain)
    end
    times = dfsums[1][:,:time]
    logz = zeros(nrow(dfsums[1]))
    ChainH5{Symbol(quant),
            typeof(quant_bma),
            typeof(times),
            typeof(logz)}(quant_bma, times, logz, 3000)
end

function bma_ss_quant(quantchains, dfsums, dirs)
    @assert length(quantchains) == length(dfsums)
    cchain = Float64[]
    weights = Float64[]
    for i in 1:length(quantchains)
        push!(cchain, quantchains[i]...)
        push!(weights, fill(exp(dfsums[i][:logz]), length(quantchains[i]))...)
    end
    return cchain, weights
end

function bmaall(dirs...)
    dfchains = Vector{DataFrame}[]
    dfsums = DataFrame[]
    models = []
    for dir in dirs
        dfchain, dfsum = load_chains(dir)
        push!(dfchains, dfchain)
        push!(dfsums, dfsum)
        push!(models, last(parsemodel(joinpath(dir, "params.dat"))))
    end
    bmamims = ROSE.StokesImage{Float64,Matrix{Float64}}[]
    for i in 1:length(dfchains[1])
        println(i)
        sums = [dfsums[j][i,:] for j in 1:length(dfsums)]
        chains = [dfchains[j][i] for j in 1:length(dfchains)]
        ims, ws = bma_ss(chains, sums, models, [dirs...])
        rims = sample(ims, Weights(ws), 500)
        push!(bmamims, mean(rims))
    end
    return bmamims
end

function bma_ss(dfchains, dfsums, models, dirs)
    @assert length(dfchains) == length(dfsums)
    @assert length(dfchains) == length(models)
    images = ROSE.StokesImage{Float64, Matrix{Float64}}[]
    weights = Float64[]
    for i in 1:length(dfchains)
        model = models[i]
        bma_images!(images, weights, dfchains[i], dfsums[i][:logz], model)
    end
    return images, weights
end

function bma_images!(images, weights, dfchain, logz, model, nimages=250)
    pr = prior(model, :amp, :cphase, :g)
    primg(x) = getfield(Soss.predict(pr, x), :img)
    for n in 1:nimages
        mimg = primg(copy(dfchain[n,3:(end-2)]))
        sim = intensitymap(mimg, 128, 128, 160.0,160.0)
        push!(images, sim)
        push!(weights, exp(logz))
    end
    return images, weights
end
