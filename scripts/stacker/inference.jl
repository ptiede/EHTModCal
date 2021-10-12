
Base.@kwdef struct RAM{T,P}
    M0::T
    n::Int
    opt_α::Float64 = 0.234
    γ::Float64 = 0.667
    q::P = Normal()
    show_progress::Bool = true
end

function optimize(l::BatchStackerLklhd, prior::NamedTuple, opt = ROSESoss.BBO())
    x0, unflatten = ParameterHandling.flatten((μ=rand(prior.μ), σ = rand(prior.σ)))
    function lpost(p)
        x = unflatten(p)
        lprior = logpdf(prior.μ, x.μ) + logpdf(prior.σ, x.σ)
        isinf(lprior) && return Inf
        return -l(x) - lprior
    end

    t = ascube(prior)
    lower = ROSESoss.transform(t, zeros(dimension(t)))
    lflat, _ = ParameterHandling.flatten(lower)

    upper = ROSESoss.transform(t, ones(dimension(t)))
    uflat, _ = ParameterHandling.flatten(upper)


    bounds = [(lflat[i], uflat[i]) for i in eachindex(lflat, uflat)]
    res = bboptimize(lpost, SearchRange=bounds, MaxFuncEvals=opt.maxevals, TraceMode=opt.tracemode)
    return unflatten(best_candidate(res)), -best_fitness(res)
end

function _getstart(p0)
    if isnothing(p0)
        x0, unflatten = ParameterHandling.flatten((μ=rand(prior.μ), σ = rand(prior.σ)))
    else
        x0, unflatten = ParameterHandling.flatten(p0)
    end
end

function StatsBase.sample(l::BatchStackerLklhd, prior::NamedTuple, s::RAM, p0=nothing)
    x0, unflatten = _getstart(p0)
    function lpost(p)
        x = unflatten(p)
        lprior = logpdf(prior.μ, x.μ) + logpdf(prior.σ, x.σ)
        isinf(lprior) && return -Inf
        return l(x) + lprior
    end
    println("Starting place ", x0)
    println("starting value ", lpost(x0))

    samples, stats, M0 = RAM_sample(lpost,
                              x0,
                              s.M0,
                              s.n;
                              opt_α=s.opt_α,
                              γ=s.γ,
                              q=s.q,
                              show_progress=s.show_progress
                            )
    tv = TupleVector([unflatten(Vector(x)) for x in eachrow(samples)])
    return tv, stats, M0
end

function StatsBase.sample(l::BatchStackerLklhd, prior::NamedTuple, s::DynestyStatic; kwargs...)
    t = ascube(prior)
    p0 = ROSESoss.transform(t, fill(0.5, dimension(t)))
    x0, unflatten = ParameterHandling.flatten(p0)
    pt(x) = first(ParameterHandling.flatten(ROSESoss.transform(t, x)))
    lp = l∘unflatten

    sampler = dynesty.NestedSampler(lp, pt, dimension(t), nlive=s.nlive)
    sampler.run_nested(;kwargs...)
    res = sampler.results
    samples, weights = res["samples"], exp.(res["logwt"] .- res["logz"][end])
    tv = TupleVector([unflatten(Vector(x)) for x in eachrow(samples)])
    return tv, weights
end

#=
Base.@kwdef struct AHMC end

function create_transform(prior)
    mins, maxs = support(prior.μ)
    TV.as((μ =TV.as(Tuple(TV.as.(Real, mins, maxs))), σ = TV.as(Tuple(TV.as.(Real, 0.0, maxs .- mins)))))
end

function StatsBase.sample(l::BatchStackerLklhd, prior::NamedTuple, s::AHMC, p0=nothing)
    t = create_transform(prior)
    function ldist(x)
        p0, ljac = TV.transform_and_logjac(t, x)
        return l(p0) + ljac
    end

end
=#
