using Distributions

function loadmodel(dir)
    model = split(dir, "/")[end]
    order = parse(Int, split(model, "-")[end])

    if occursin("mring_floor_order", model)
        img = ROSESoss.mringwfloor(N=order,)
    elseif occursin("mring_gfloor_order", model)
        img = ROSESoss.mringwgfloor(N=order,)
    elseif occursin("mring_gfloor_newprior_order", model)
        img = ROSESoss.mringwgfloor(N=order,)
    elseif occursin("mring_order", model)
        img = ROSESoss.mring(N=order,)
    elseif occursin("mring_floor_stretch_order", model)
        img = ROSESoss.smringwfloor(N=order,)
    else
        throw("model $model not found")
    end
    return img
end

function rand!(f, samples, dist)
    for I in CartesianIndices(samples)
        samples[I] = f(rand(dist))
    end
end


function transform_sparam(f::T, means, stds, mins, maxs; nsamples=25_000) where {T}
    tmean = zeros(length(means))
    tstd = zeros(length(stds))
    samples = zeros(nsamples)
    for i in eachindex(means, stds)
        dist = product_distribution(truncated.(Normal.(means[i], stds[i]), mins, maxs))
        rand!(f, samples, dist)
        tmean[i], tstd[i] = mean_and_std(samples)
    end
    return tmean, tstd
end

function construct_meanstd(df, labels::Vector{String})
    means = [Tuple([df[i, "μ_"*l] for l in labels]) for i in 1:nrow(df)]
    std = [Tuple([df[i, "σ_"*l] for l in labels]) for i in 1:nrow(df)]
    return means, std
end

function floorfluxfrac(x::AbstractVector)
    return @inbounds floorfluxfrac(x[1], x[2], x[3], x[4], x[5])
end

function floorfluxfrac(flux, diam, fwhm, floor, dg)
    σg = dg/(2*sqrt(2*log(2)))
    r = diam/2
    fg = floor*flux
    ff = fg*(1 - exp(-(r/σg)^2/2))
    return ff/(ff + (1-floor)*flux)
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

function convertall(dir)
    sdirs = filter(isdir, readdir(dir, join=true))
    cfiles = joinpath.(sdirs, Ref("chain.h5"))
    convert_antonio.(cfiles)
end

function convert_antonio(cfile, tmin=12.455, tmax=14.2)
    println("Converting $cfile")
    chain = ChainH5(cfile, (:img_diam, :img_fwhm, :img_mp_1))
    btchain = restricttime(chain, tmin, tmax)
    save_antonio(joinpath(dirname(cfile), "chain.txt"), btchain)
end

function save_antonio(out, chain::ChainH5)
    diam = getparam(chain, :img_diam).chain
    α = getparam(chain, :img_fwhm).chain
    pa = getparam(chain, :img_mp_1).chain
    name,ext = splitext(out)
    times = chain.times
    open(name*"_diam"*ext, "w") do io
        for i in eachindex(times, diam)
            pdiam = @. diam[i] - 1/(8*log(2))*α[i]^2/diam[i]
            writedlm(io, [times[i] pdiam'])
        end
    end

    open(name*"_pa"*ext, "w") do io
        for i in eachindex(times, diam)
            tmp = -rad2deg.(pa[i])
            writedlm(io, [times[i] tmp'])
        end
    end
end


function parsechainpath(cfile)
    model = split(cfile, "/")[end]
    mins = [25.0, #diam
            1.0, #width
            0.1, #flux
            ]
    maxs = [85.0,
            40.0,
            4.0
           ]
    quants = [:img_diam, :img_fwhm, :img_f]
    labels = ["diameter (μas)", "width (μas)", "flux (Jy)"]
    wrapped = Bool[false, false, false]

    if occursin("gfloor", model)
        push!(mins, 0.0, 40.0)
        push!(maxs, 1.0, 200.0)
        push!(quants, :img_floor, :img_dg)
        push!(wrapped, false, false)
        push!(labels, "floor flux fraction", "Gaussian diameter")
    elseif occursin("floor", model)
        push!(mins, 0.0)
        push!(maxs, 1.0)
        push!(quants, :img_floor)
        push!(labels, "floor flux fraction")
        push!(wrapped, false)
    end

    if occursin("stretch", model)
        push!(mins, 0.0, -π/2)
        push!(maxs, 0.5, π/2)
        push!(quants, :img_τ, :img_ξτ)
        push!(labels, "ellipticity τ", "ellipticity PA (rad) W of N")
        push!(wrapped, false, true)
    end

    order = parse(Int, split(model, "-")[end])

    for i in 1:order
        push!(quants, Symbol("img_ma_$i"), Symbol("img_mp_$i"))
        push!(mins, 0.0, -1π)
        push!(maxs, 0.5, 1π)
        push!(labels, "m=$i amp", "m=$i phase (deg) E of N")
        push!(wrapped, false, true)
    end

    return mins, maxs, wrapped, Tuple(quants), labels

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


function df2tv(df::DataFrame)
    nms = filter(x->!(occursin("σ_img_", x)), names(df))
    iimg = findall(x->occursin("μ_img_",x), nms)
    inon = findall(x->!((occursin("μ_img_mp_",x)||occursin("μ_img_ma_",x))&&occursin("μ_img_", x)), nms)
    order = length(findall(x->occursin("μ_img_ma_",x), nms))
    ks = Symbol.(last.(split.(nms[inon], Ref("μ_img_"))))
    ntproto = namedtuple(ks..., :ma, :mp)
    ma = zeros(order, nrow(df))
    mp = zeros(order, nrow(df))
    for i in 1:order
        ma[i,:] = getproperty(df, Symbol("μ_img_ma_$i"))
        mp[i,:] = getproperty(df, Symbol("μ_img_mp_$i"))
    end
    imgt = Tuple(getproperty(df, Symbol(k)) for k in nms[inon])
    return TupleVector(ntproto(imgt..., nestedview(ma), nestedview(mp)))
end
