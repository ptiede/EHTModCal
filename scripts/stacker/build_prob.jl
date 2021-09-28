function build_model_all(θ::NamedTuple, mins, maxs, wrapped)
    μ, σ = θ.μ, θ.σ
    dists = Union{Normal{Float64}, EHTModelStacker.VonMisesWrap{Float64, Float64}}[]
    for i in eachindex(wrapped)
        if wrapped[i]
            push!(dists, EHTModelStacker.VonMisesWrap(μ[i], σ[i]))
        else
            push!(dists, Normal(μ[i], σ[i]))
        end
    end
    transition = EHTModelStacker.MyProduct(dists)#EHTModelStacker.MvNormalFast(μ, Σ.^2)
    prior = EHTModelStacker.MvUniform(mins, maxs)
    return SnapshotWeights(transition, prior)
end


function build_loglklhd(build_model, chain, min, max, wrapped)
    return function (θ::NamedTuple)
        ws = build_model(θ, min, max, wrapped)
        return lpdf(ws, chain)
    end
end
