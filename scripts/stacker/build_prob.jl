function EHTModelStacker.SnapshotWeights(θ, mins, maxs, wrapped, batchsize)
    μ, σ = θ.μ, θ.σ
    dists = Union{Truncated{EHTModelStacker.NormalFast{Float64}, Continuous, Float64}, EHTModelStacker.VonMisesWrap{Float64, Float64}}[]
    for i in eachindex(wrapped)
        if wrapped[i]
            push!(dists, EHTModelStacker.VonMisesWrap(μ[i], σ[i]))
        else
            push!(dists, truncated(EHTModelStacker.NormalFast(μ[i], σ[i]), mins[i], maxs[i]))
        end
    end
    #dists = EHTModelStacker.NormalFast.(μ, σ)
    transition = EHTModelStacker.MyProduct(dists)#EHTModelStacker.MvNormalFast(μ, Σ.^2)
    prior = EHTModelStacker.MvUniform(mins, maxs)
    return SnapshotWeights(transition, prior, batchsize)
end


#=
function EHTModelStacker.SnapshotWeights(θ, mins, maxs, wrapped, batchsize)
    μ, σ = θ.μ, θ.σ
    dists = Union{EHTModelStacker.NormalFast{Float64}, EHTModelStacker.VonMisesWrap{Float64, Float64}}[]
    for i in eachindex(wrapped)
        if wrapped[i]
            push!(dists, EHTModelStacker.VonMisesWrap(μ[i], σ[i]))
        else
            push!(dists, EHTModelStacker.NormalFast(μ[i], σ[i]))
        end
    end
    #dists = EHTModelStacker.NormalFast.(μ, σ)
    transition = EHTModelStacker.MyProduct(dists)#EHTModelStacker.MvNormalFast(μ, Σ.^2)
    prior = EHTModelStacker.MvUniform(mins, maxs)
    return SnapshotWeights(transition, prior, batchsize)
end
=#

struct BatchStackerLklhd{C, T, B}
    chain::C
    min::T
    max::T
    wrapped::B
    batchsize::Int
end

function (l::BatchStackerLklhd)(θ)
    ws = SnapshotWeights(θ, l.min, l.max, l.wrapped, l.batchsize)
    lapprox = lpdf(ws, l.chain)
    return lapprox
end
