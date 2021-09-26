function average_chain_ram(build_model, minp, maxp, chain; nsamples=10_000)
    pμ = Uniform(minp, maxp)
    pσ = Uniform(0.0, 20.0)
    p0 = (μ=rand(pμ), σ=rand(pσ))
    lprior(x) = logpdf(pμ, x[1]) + logpdf(pσ, x[2])
    fp0, unflatten = ParameterHandling.flatten(p0)
    println(fp0)
    llklhd = build_loglklhd(build_model, chain, minp, maxp)∘unflatten
    function lpost(x)
        lp = lprior(x)
        isinf(lp) ? -1e300 : llklhd(x) + lp
    end
    chain, acc = RAM_sample(lpost, [50.0, 0.2], [0.1,0.2], nsamples; show_progress=true)
    return TupleVector((μ=chain[nsamples÷2:end,1], σ=chain[nsamples÷2:end,2]))
end


function build_model_LN(θ::NamedTuple, min, max)
    μ, σ = θ.μ, θ.σ
    transition = LogNormal(log(μ), σ)
    prior = Uniform(min, max)
    return SnapshotWeights(transition, prior)
end

function build_model_N(θ::NamedTuple, min, max)
    μ, σ = θ.μ, θ.σ
    transition = Normal(μ, σ)
    prior = Uniform(min, max)
    return SnapshotWeights(transition, prior)
end

function build_model_C(θ::NamedTuple, min, max)
    μ, σ = θ.μ, θ.σ
    transition = Cauchy(μ, σ)
    prior = Uniform(min, max)
    return SnapshotWeights(transition, prior)
end

function build_model_K(θ::NamedTuple, min, max, N)
    knots = range(min, max, length=N)
    w, bw = θ.w, θ.bw
    transition = kde(knots, bandwidth=bw, boundary=(min/2, max*2), weights=Weights(w))
    prior = Uniform(min, max)
    return SnapshotWeights(transition, prior)
end



function build_model_dw(θ::NamedTuple, min, max)
    μ, σ, ρ = θ.μ, θ.σ, θ.ρ
    s2 = σ[1]^2
    t2 = σ[2]^2
    c = σ[1]*σ[2]*ρ
    cov = @SMatrix [s2 c; c t2]
    transition = MvNormal2D(μ, cov)
    prior = product_distribution([Uniform(min[1],max[1]), Uniform(min[2],max[2])])
    return SnapshotWeights(transition, prior)
end




function average_chain(build_model, minp, maxp, chain)
    p0 = (μ=(maxp-minp/2)+minp, σ=(maxp-minp)/2)
    fp0, unflatten = ParameterHandling.flatten(p0)
    lp = build_loglklhd(build_model, chain, minp, maxp)∘unflatten
    pt = build_prior_transform(minp, maxp)

    fu = Dict("min_ncall"=>100, "min_eff"=>0.2)
    sampler = dynesty.NestedSampler(lp, pt, 2, first_update=fu)
    sampler.run_nested()
    res = sampler.results
    samples, weights = res["samples"], exp.(res["logwt"] .- res["logz"][end])
    tv = TupleVector((μ=samples[:,1], σ=samples[:,2]))
    echain = TupleVector(sample(tv, Weights(weights), 2000))
    return echain
end

function average_chain_dw(build_model, minp, maxp, chain)
    priors = (μ = Product(Uniform.(minp,maxp)), σ=Product(Uniform.(0.1,(maxp.-minp))), ρ=Uniform(0.0, 0.99))
    t = ascube(priors)
    p0 = HypercubeTransform.transform(t, fill(0.5, dimension(t)))
    fp0, unflatten = ParameterHandling.flatten(p0)
    lp = build_loglklhd(build_model, chain, minp, maxp)∘unflatten
    pt(x) = first(ParameterHandling.flatten(HypercubeTransform.transform(t, x)))
    #return lp, pt, unflatten
    sampler = dynesty.NestedSampler(lp, pt, dimension(t))
    sampler.run_nested(dlogz=1.0)
    res = sampler.results
    samples, weights = res["samples"], exp.(res["logwt"] .- res["logz"][end])
    tv = TupleVector([unflatten(Vector(x)) for x in eachrow(samples)])
    echain = TupleVector(sample(tv, Weights(weights), 2000))
    return echain
end
