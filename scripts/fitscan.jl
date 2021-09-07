function _run_model(logj, jscan, sampler; kwargs...)
    chain, stats = sampler(logj; kwargs...)
    logz = stats[:logz]
    logzerr = stats[:logzerr]
    logl = logl

    sess = 1/sum(abs2, chain.weights)
    echain = TupleVector(sample(chain, Weights(vec(chain.weights)), 2*Int(floor(sess))))
    @info "Nested sampling finished ess = $(sess)"

    DataFrame(chain) |> CSV.write(joinpath(@__DIR__,dir, "Chain",  "nested_chain_scan-$jscan.csv"))
    DataFrame(echain) |> CSV.write(joinpath(@__DIR__,dir, "Chain",  "equal_chain_scan-$jscan.csv"))

    # Construct the image
    opt = delete(NamedTuple{pnames}(Array(chain[end])), :weights)
    mopt = Soss.predict(logj.model, opt)[:img]
    return logz, logzerr, logl, opt, mopt
end


function fit_scan(dir::String,
                              model,
                              ampobs::ROSE.EHTObservation{F,A},
                              cpobs::ROSE.EHTObservation{F,P},
                              jscan::Int,
                              sampler;
                              kwargs...
                              ) where {F,A<:ROSE.EHTVisibilityAmplitudeDatum,
                                         P<:ROSE.EHTClosurePhaseDatum}

    tstart = round(ampobs[1].time, digits=2)

    logj = ROSESoss.create_joint(model, ampobs, cpobs)
    tc = ascube(logj)
    ndim = dimension(tc)

    bl = getdata(ampobs,:bl)
    s1 = unique(first.(bl))
    s2 = unique(last.(bl))
    stations = unique([s1...,s2...])


    logz, logzerr, logl, opt, mopt = _run_model(logj, jscan, sampler; kwargs...)


    gs = Soss.predict(logj.model, opt)[:g]
    achi2 = ROSESoss.chi2(mopt, ampobs, gs)
    cpchi2 = ROSESoss.chi2(mopt, cpobs)
    rchi2 = (achi2+cpchi2)/(ROSE.nsamples(ampobs)+ROSE.nsamples(cpobs)-ndim)
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-ndim)
    rcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs)-(ndim-length(stations)))
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    mcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs))
    @info "Reduced chi square: $rchi2"

    return merge((time = tstart, rchi2=rchi2, rampchi2=rachi2, rcpchi2=rcpchi2,
                    mampchi2=machi2, mcpchi2=mcpchi2,
                    logz=logz, logzerr=logzerr, logl=logl), opt)
end

function fit_scan(dir::String,
    model,
    ampobs::ROSE.EHTObservation{F,A},
    cpobs::ROSE.EHTObservation{F,P},
    jscan::Int,
    sampler;
    fitgains=false,
    plot_results=true,
    obschar=ROSESoss.ObsChar(),
    kwargs...
    ) where {F,A<:ROSE.EHTLogClosure,
               P<:ROSE.EHTClosurePhaseDatum}

    tstart = round(ampobs[1].time, digits=2)

    logj = ROSESoss.create_joint(model, ampobs, cpobs)
    tc = ascube(logj)
    ndim = dimension(tc)

    logz, logzerr, logl, opt, mopt = _run_model(logj, jscan, sampler; kwargs...)

    achi2 = ROSESoss.chi2(mopt, ampobs)
    cpchi2 = ROSESoss.chi2(mopt, cpobs)
    rchi2 = (achi2+cpchi2)/(ROSE.nsamples(ampobs)+ROSE.nsamples(cpobs)-ndim)
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-ndim)
    rcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs)-ndim)
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    mcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs))
    @info "Reduced chi square: $rchi2"


    return merge((time = tstart, rchi2=rchi2, rampchi2=rachi2, rcpchi2=rcpchi2,
                    mampchi2=machi2, mcpchi2=mcpchi2,
                    logz=logz, logzerr=logzerr, logl=logl), opt)
end

function fit_scan(dir::String,
                              model,
                              ampobs::ROSE.EHTObservation{F,A},
                              jscan::Int,
                              sampler;
                              kwargs...
                              ) where {F,A<:ROSE.EHTVisibilityAmplitudeDatum}


    tstart = round(ampobs[1].time, digits=2)
    logj = ROSESoss.create_joint(model, ampobs)
    tc = ascube(logj)
    ndim = dimension(tc)


    logz, logzerr, logl, opt, mopt = _run_model(logj, jscan, sampler; kwargs...)

    achi2 = ROSESoss.chi2(mopt, ampobs, gs)
    rchi2 = (achi2)/(ROSE.nsamples(ampobs)-ndim)
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-ndim)
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    @info "Reduced chi square: $rchi2"

    return merge((time = tstart, rchi2=rchi2, rampchi2=rachi2,
                    mampchi2=machi2,
                    logz=logz, logzerr=logzerr, logl=logl), opt)
end

function fit_scan(dir::String,
                              model,
                              visobs::ROSE.EHTObservation{F,V},
                              jscan::Int, sampler;
                              kwargs...
                              ) where {F,V<:ROSE.EHTVisibilityDatum}

    tstart = round(ampobs[1].time, digits=2)
    logj = ROSESoss.create_joint(model, visobs)
    tc = ascube(logj)
    ndim = dimension(tc)

    logz, logzerr, logl, opt, mopt = _run_model(logj, jscan, sampler; kwargs...)



    vchi2 = ROSESoss.chi2(mopt, visobs, ga, gp)
    rchi2 = vchi2/(2*ROSE.nsamples(visobs)-ndim)
    mchi2 = (vchi2)/(2*ROSE.nsamples(visobs))
    @info "Reduced chi square: $rchi2"


    return merge((time = tstart, rchi2=rchi2,
                    mchi2=mchi2,
                    logz=logz, logzerr=logzerr, logl=logl), opt)
end
