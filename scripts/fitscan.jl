function _run_model(dir, logj, sampler; kwargs...)
    chain, stats = sampler(logj; kwargs...)
    logz = stats[:logz]
    logzerr = stats[:logzerr]
    logl = stats[:logl]

    sess = 1/sum(abs2, chain.weights)
    echain = TupleVector(sample(chain, Weights(vec(chain.weights)), 2*Int(floor(sess))))
    @info "Nested sampling finished ess = $(sess)"

    # Construct the image
    opt = delete(chain[end], :weights)
    model = Soss.argvals(logj)[:image]
    mopt = Soss.predict(model, opt[:img])
    return logz, logzerr, logl, opt, mopt, chain, echain
end

function nt2df(nt::NamedTuple)
    df = DataFrame()
    _nt2df!(df, nt)
    df
end

function _nt2df!(df, nt::NamedTuple)
    for k in keys(nt)
        _nt2df!(df, nt[k], String(k))
    end
end

function _nt2df!(df, nt::NamedTuple, name::String)
    for k in keys(nt)
        tmp = name*"_"*String(k)
        _nt2df!(df, nt[k], tmp)
    end
end

function _nt2df!(df, nt::Number, name::String)
    insertcols!(df, name=>nt)
end

function _nt2df!(df, nt::AbstractVector, name::String)
    for (i,v) in enumerate(nt)
        insertcols!(df, Symbol(name*"_"*"$i") => v)
    end
end

function tv2df(tv::TupleVector)
    df = DataFrame()
    _tv2df!(df, tv)
    df
end

function _tv2df!(df, tv::TupleVector)
    for k in propertynames(tv)
        name = String(k)
        _tv2df!(df, getproperty(tv, Symbol(k)), name)
    end
end

function _tv2df!(df, tv::TupleVector, name::String)
    for k in propertynames(tv)
        tmp = name*"_"*String(k)
        _tv2df!(df, getproperty(tv, Symbol(k)), tmp)
    end
end

function _tv2df!(df, tv::AbstractVector, name::String)
    if typeof(tv[1]) <: Number
        insertcols!(df, Symbol(name) => tv)
    else
        s = length(tv[1])
        for i in 1:s
            insertcols!(df, Symbol(name*"_"*"$i") => getindex.(tv,i))
        end
    end
end



function fit_scan(dir::String,
                              model,
                              ampobs::ROSE.EHTObservation{F,A},
                              cpobs::ROSE.EHTObservation{F,P},
                              sampler;
                              kwargs...
                              ) where {F,A<:ROSE.EHTVisibilityAmplitudeDatum,
                                         P<:ROSE.EHTClosurePhaseDatum}

    tstart = round(ampobs[1].time, digits=2)

    logj = ROSESoss.create_joint(model, ampobs, cpobs)
    tc = ascube(logj)
    ndim = dimension(tc)

    bl = ROSE.getdata(ampobs,:baselines)
    s1 = unique(first.(bl))
    s2 = unique(last.(bl))
    stations = unique([s1...,s2...])


    logz, logzerr, logl, opt, mopt, chain, echain = _run_model(dir, logj, sampler; kwargs...)


    gs = Soss.predict(argvals(logj)[:gamps], opt[:g])
    achi2 = ROSESoss.chi2(mopt, ampobs, gs)
    cpchi2 = ROSESoss.chi2(mopt, cpobs)
    rchi2 = (achi2+cpchi2)/(ROSE.nsamples(ampobs)+ROSE.nsamples(cpobs)-ndim)
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-ndim)
    rcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs)-(ndim-length(stations)))
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    mcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs))
    @info "Reduced chi square: $rchi2"
    stats = merge((time = tstart, rchi2=rchi2, rampchi2=rachi2, rcpchi2=rcpchi2,
                   mampchi2=machi2, mcpchi2=mcpchi2,
                   logz=logz, logzerr=logzerr, logl=logl), opt)


    return chain, echain, stats
end




function fit_scan(dir::String,
    model,
    ampobs::ROSE.EHTObservation{F,A},
    cpobs::ROSE.EHTObservation{F,P},
    sampler;
    kwargs...
    ) where {F,A<:ROSE.EHTLogClosureAmplitudeDatum,
               P<:ROSE.EHTClosurePhaseDatum}

    tstart = round(ampobs[1].time, digits=2)

    logj = ROSESoss.create_joint(model, ampobs, cpobs)
    tc = ascube(logj)
    ndim = dimension(tc)

    logz, logzerr, logl, opt, mopt, chain, echain = _run_model(dir, logj, sampler; kwargs...)

    achi2 = ROSESoss.chi2(mopt, ampobs)
    cpchi2 = ROSESoss.chi2(mopt, cpobs)
    rchi2 = (achi2+cpchi2)/(ROSE.nsamples(ampobs)+ROSE.nsamples(cpobs)-ndim)
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-ndim)
    rcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs)-ndim)
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    mcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs))
    @info "Reduced chi square: $rchi2"


    stats = merge((time = tstart, rchi2=rchi2, rampchi2=rachi2, rcpchi2=rcpchi2,
                    mampchi2=machi2, mcpchi2=mcpchi2,
                    logz=logz, logzerr=logzerr, logl=logl), opt)
    return chain, echain, stats
end

function fit_scan(dir::String,
                              model,
                              ampobs::ROSE.EHTObservation{F,A},
                              sampler;
                              kwargs...
                              ) where {F,A<:ROSE.EHTVisibilityAmplitudeDatum}


    tstart = round(ampobs[1].time, digits=2)
    logj = ROSESoss.create_joint(model, ampobs)
    tc = ascube(logj)
    ndim = dimension(tc)


    logz, logzerr, logl, opt, mopt, chain, echain = _run_model(dir, logj, sampler; kwargs...)

    achi2 = ROSESoss.chi2(mopt, ampobs, gs)
    rchi2 = (achi2)/(ROSE.nsamples(ampobs)-ndim)
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-ndim)
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    @info "Reduced chi square: $rchi2"

    stats = merge((time = tstart, rchi2=rchi2, rampchi2=rachi2,
                    mampchi2=machi2,
                    logz=logz, logzerr=logzerr, logl=logl), opt)
    return chain, echain, stats

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

    logz, logzerr, logl, opt, mopt, chain, echain = _run_model(dir, logj, sampler; kwargs...)



    vchi2 = ROSESoss.chi2(mopt, visobs, ga, gp)
    rchi2 = vchi2/(2*ROSE.nsamples(visobs)-ndim)
    mchi2 = (vchi2)/(2*ROSE.nsamples(visobs))
    @info "Reduced chi square: $rchi2"


    stats = merge((time = tstart, rchi2=rchi2,
                    mchi2=mchi2,
                    logz=logz, logzerr=logzerr, logl=logl), opt)

    return chain, echain, stats

end
