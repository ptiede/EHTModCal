using HDF5
using DataFrames, CSV
using StatsBase

function getparamnames(fname)
    tmp = split(fname, "_")
    day = tmp[2]
    mflux = tmp[3]
    spin = parse(Float64, tmp[4][2:end])
    rhigh = parse(Float64, tmp[5][3:end])
    inclination = parse(Float64, tmp[6][2:end])
    band = tmp[end]
    return (;day, band, mflux, rhigh, spin, inclination)
end


function createh5(flist, truths, order, outname, nsamples=1000)
    samples = zeros(length(flist), nsamples)
    for (i,f) in enumerate(flist)
        df = CSV.read(joinpath(f, "snapshot_120/noisefrac0.02/mring_gfloor_order-$order", "ChainHA_3/chain_ha_trunc.csv"), DataFrame)
        s = df[end÷2:end, :μ_img_diam]
        samples[i, :] .= sample(s, nsamples)
    end
    pnames = DataFrame(getparamnames.(flist))
    addscale!(pnames, truths)
    println(pnames)
    h5open(outname, "w") do fid
        for n in names(pnames)
            write(fid, n, pnames[:,n])
        end
        write(fid, "samples", samples)
    end
    return samples
end


function addscale!(pnames, tfile)
    truths = CSV.read(tfile, DataFrame)
    pnames.thetag = zeros(nrow(pnames))
    for i in 1:nrow(pnames)
        for j in 1:nrow(truths)
            if ((pnames[i,:spin] == truths[j,:spin])&&
                (pnames[i,:inclination]==truths[j,:inc])&&
                (pnames[i,:rhigh]==truths[j,:Rhi])&&
                (pnames[i,:mflux]==truths[j,:acc]))
                pnames[i,:thetag] = truths[j,:scale]*5.03
            end
        end
    end
end
