
Base.@kwdef struct RAM{T,P}
    M0::T = 0.1
    n::Int = 100_000
    opt_α::Float64 = 0.234
    γ::Float64 = 0.667
    q::P = Normal()
    show_progress::Bool = true
end

function ROSESoss.optimize(l::BatchStackerLklhd, prior::NamedTuple, opt = ROSESoss.BBO(;maxevals=20_000))
    t = ascube(prior)
    lp(x) = -l(HypercubeTransform.transform(t, x))
    lp(rand(dimension(t)))

    bounds = [(0.0, 1.0) for i in 1:dimension(t)]
    res = bboptimize(lp, SearchRange=bounds, MaxFuncEvals=opt.maxevals, TraceMode=opt.tracemode)
    return best_candidate(res), -best_fitness(res)
end

function _getstart(p0, prior)
    if isnothing(p0)
        x0, unflatten = ParameterHandling.flatten((μ=rand(prior.μ), σ = rand(prior.σ)))
    else
        x0, unflatten = ParameterHandling.flatten(p0)
    end
end

function StatsBase.sample(l::BatchStackerLklhd, pr::NamedTuple, s::RAM, p0)
    t = ascube(pr)

    # The famed julia closure bug!
    function lp(x)
        for i in x
            if !(0.0 < i < 1.0)
                return -Inf
            end
        end
        l(HypercubeTransform.transform(t, x))
    end
    #compile the function because of julia bug I guess?
    println("Starting place ", p0)
    println("starting value ", lp(p0))
    samples, stats, M0,_ = RAM_sample(lp,
                              rand(dimension(t)),
                              s.M0,
                              s.n;
                              opt_α=s.opt_α,
                              γ=s.γ,
                              q=s.q,
                              show_progress=s.show_progress
                            )
    tv = TupleVector([HypercubeTransform.transform(t, x) for x in eachrow(samples)])
    return tv, stats, M0
end

function StatsBase.sample(l::BatchStackerLklhd, prior::NamedTuple, s::DynestyStatic; kwargs...)
    t = ascube(prior)
    p0 = ROSESoss.transform(t, fill(0.5, dimension(t)))
    x0, unflatten = ParameterHandling.flatten(p0)
    pt(x) = first(ParameterHandling.flatten(ROSESoss.transform(t, x)))
    lp = l∘unflatten

    sampler = ROSESoss.dynesty.NestedSampler(lp, pt, dimension(t), nlive=s.nlive)
    sampler.run_nested(;kwargs...)
    res = sampler.results
    samples, weights = res["samples"], exp.(res["logwt"] .- res["logz"][end])
    tv = TupleVector([unflatten(Vector(x)) for x in eachrow(samples)])
    return sample(tv, Weights(weights), 10_000)
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
