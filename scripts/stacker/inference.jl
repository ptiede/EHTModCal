
function ROSESoss.optimize(l::BatchStackerLklhd, prior::NamedTuple, opt = ROSESoss.BBO(;maxevals=40_000))
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

function StatsBase.sample(l::BatchStackerLklhd, pr::NamedTuple, s::RAM, nsteps; kwargs...)
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
    #compile the function because of julia bug I guess
    p0 = s.x
    println("Starting place ", p0)
    println("starting value ", lp(p0))
    samples, stats, state,logprob = sample(lp,
                                           s, nsteps;
                                           kwargs...
                                          )
    tv = TupleVector([HypercubeTransform.transform(t, x) for x in eachrow(samples)])
    return tv, stats, state, logprob
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
