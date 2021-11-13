using Pkg; Pkg.activate(@__DIR__)
using HDF5
using DataFrames
using CSV


function make_h5_summary(files, hfile, parameter, truths)
    h5open(file, "W") do fid
        create_group
    end
end
