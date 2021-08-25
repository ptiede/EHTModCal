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
                              ) where {F,A<:ROSE.EHTVisibilityAmplitudeDatum,
                                         P<:ROSE.EHTClosurePhaseDatum}

    tstart = round(ampobs[1].time, digits=2)
    bl = getdata(ampobs,:baselines)
    s1 = unique(first.(bl))
    s2 = unique(last.(bl))
    stations = unique([s1...,s2...])

    logj = ROSESoss.create_joint(model, ampobs, cpobs; fitgains=fitgains)

    chain, state, pnames = sampler(logj; kwargs...)
    logz = chain.logevidence
    logzerr = state.logzerr

    sess = 1/sum(abs2, chain[:weights])
    echain = sample(chain, Weights(vec(chain[:weights])), 2*Int(floor(sess)))
    @info "Nested sampling finished ess = $(sess)"

    DataFrame(echain) |> CSV.write(joinpath(@__DIR__,dir, "Chain",  "nested_chain_scan-$jscan.csv"))

    # Construct the image
    opt = NamedTuple{pnames}(Array(chain[end]))
    mopt = Soss.predict(logj.model, opt)[:img]
    #println("here1")
    #intensitymap(mopt, 120,120,120,120)
    #println("here2")

    #Find a chi-square because yes...
    if fitgains
        gs = (AP=opt[:AP], AZ=opt[:AZ], JC=opt[:JC], SM=opt[:SM], AA=opt[:AA],
                  LM=opt[:LM], SP=opt[:SP], PV=opt[:PV])
    else
        gs = (AP=1.0, AZ=1.0, JC=1.0, SM=1.0, AA=1.0,
                  LM=1.0, SP=1.0, PV=1.0)
    end
    achi2 = ROSESoss.chi2(mopt, ampobs, gs)
    cpchi2 = ROSESoss.chi2(mopt, cpobs)
    rchi2 = (achi2+cpchi2)/(ROSE.nsamples(ampobs)+ROSE.nsamples(cpobs)-(length(opt)-(length(gs)-length(stations))))
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-(length(opt)-(length(gs)-length(stations))))
    rcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs)-length(opt)+length(gs))
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    mcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs))
    @info "Reduced chi square: $rchi2"

    #write summary files
    if plot_results
        samples = Array(echain)

        params = [NamedTuple{pnames}(samples[i,:]) for i in rand(1:size(samples,1), 500)]
        ims = ROSESoss.make_mean(logj.model, params)
        p = ROSESoss.plot_im(ims)
        title!(p, "Time: $(tstart) UTC")
        savefig(joinpath(@__DIR__, dir,  "MAP","mean_fit_$jscan.png"))

        #Plot some fit statistics stuff for MAP
        p1=plot_amp_comp(ampobs, gs, mopt)
        title!(p1,"Time: $tstart   χ²ᵣ = $rachi2")
        savefig(p1,joinpath(@__DIR__, dir,  "Residual", "amp_residuals_scan-$jscan.png"))

        p2=plot_cp_comp(cpobs, mopt)
        title!(p2,"Time: $tstart))   χ²ᵣ = $rcpchi2")
        savefig(p2, joinpath(@__DIR__, dir,  "Residual", "cp_residuals_scan-$jscan.png"))
        Plots.closeall()
    end

    return merge((time = tstart, rchi2=rchi2, rampchi2=rachi2, rcpchi2=rcpchi2,
                    mampchi2=machi2, mcpchi2=mcpchi2,
                    logz=logz, logzerr=logzerr), opt)
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
    bl = getdata(ampobs,:baselines)
    s1 = unique(first.(bl))
    s2 = unique(last.(bl))
    stations = unique([s1...,s2...])

    logj = ROSESoss.create_joint(model, ampobs, cpobs; fitgains=fitgains)

    chain, state, pnames = sampler(logj; kwargs...)
    logz = chain.logevidence
    logzerr = state.logzerr

    sess = 1/sum(abs2, chain[:weights])
    echain = sample(chain, Weights(vec(chain[:weights])), 2*Int(floor(sess)))
    @info "Nested sampling finished ess = $(sess)"

    DataFrame(echain) |> CSV.write(joinpath(@__DIR__,dir, "Chain",  "nested_chain_scan-$jscan.csv"))

    # Construct the image
    opt = NamedTuple{pnames}(Array(chain[end]))
    mopt = Soss.predict(logj.model, opt)[:img]
    #println("here1")
    #intensitymap(mopt, 120,120,120,120)
    #println("here2")

    #Find a chi-square because yes...
    if fitgains
    gs = (AP=opt[:AP], AZ=opt[:AZ], JC=opt[:JC], SM=opt[:SM], AA=opt[:AA],
    LM=opt[:LM], SP=opt[:SP], PV=opt[:PV])
    else
    gs = (AP=1.0, AZ=1.0, JC=1.0, SM=1.0, AA=1.0,
    LM=1.0, SP=1.0, PV=1.0)
    end
    achi2 = ROSESoss.chi2(mopt, ampobs, gs)
    cpchi2 = ROSESoss.chi2(mopt, cpobs)
    rchi2 = (achi2+cpchi2)/(ROSE.nsamples(ampobs)+ROSE.nsamples(cpobs)-(length(opt)-(length(gs)-length(stations))))
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-(length(opt)-(length(gs)-length(stations))))
    rcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs)-length(opt)+length(gs))
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    mcpchi2 = (cpchi2)/(ROSE.nsamples(cpobs))
    @info "Reduced chi square: $rchi2"

    #write summary files
    if plot_results
    samples = Array(echain)

    params = [NamedTuple{pnames}(samples[i,:]) for i in rand(1:size(samples,1), 500)]
    ims = ROSESoss.make_mean(logj.model, params)
    p = ROSESoss.plot_im(ims)
    title!(p, "Time: $(tstart) UTC")
    savefig(joinpath(@__DIR__, dir,  "MAP","mean_fit_$jscan.png"))

    #Plot some fit statistics stuff for MAP
    p1=plot_amp_comp(ampobs, gs, mopt)
    title!(p1,"Time: $tstart   χ²ᵣ = $rachi2")
    savefig(p1,joinpath(@__DIR__, dir,  "Residual", "amp_residuals_scan-$jscan.png"))

    p2=plot_cp_comp(cpobs, mopt)
    title!(p2,"Time: $tstart))   χ²ᵣ = $rcpchi2")
    savefig(p2, joinpath(@__DIR__, dir,  "Residual", "cp_residuals_scan-$jscan.png"))
    Plots.closeall()
    end

    return merge((time = tstart, rchi2=rchi2, rampchi2=rachi2, rcpchi2=rcpchi2,
    mampchi2=machi2, mcpchi2=mcpchi2,
    logz=logz, logzerr=logzerr), opt)
end

function fit_scan(dir::String,
                              model,
                              ampobs::ROSE.EHTObservation{F,A},
                              jscan::Int,
                              sampler;
                              fitgains=false,
                              plot_results=true,
                              obschar=ROSESoss.ObsChar(),
                              kwargs...
                              ) where {F,A<:ROSE.EHTVisibilityAmplitudeDatum}

    tstart = round(ampobs[1].time, digits=2)
    bl = getdata(ampobs,:bl)
    s1 = unique(first.(bl))
    s2 = unique(last.(bl))
    stations = unique([s1...,s2...])

    logj = ROSESoss.create_joint(model, ampobs; fitgains=fitgains)

    chain, state, pnames = sampler(logj; kwargs...)
    logz = chain.logevidence
    logzerr = state.logzerr

    ess = 1/sum(abs2, chain[:weights])
    echain = sample(chain, Weights(vec(chain[:weights])), 2*Int(floor(ess)))
    @info "Nested sampling finished ess = $(ess)"

    DataFrame(echain) |> CSV.write(joinpath(@__DIR__,dir, "Chain",  "nested_chain_scan-$jscan.csv"))

    # Construct the image
    opt = NamedTuple{pnames}(Array(chain[end]))
    mopt = Soss.predict(logj.model, opt)[:img]
    intensitymap(mopt, 120,120,120,120)

    #Find a chi-square because yes...
    if fitgains
        gs = (AP=opt[:AP], AZ=opt[:AZ], JC=opt[:JC], SM=opt[:SM], AA=opt[:AA],
                  LM=opt[:LM], SP=opt[:SP])
    else
        gs = (AP=1.0, AZ=1.0, JC=1.0, SM=1.0, AA=1.0,
                  LM=1.0, SP=1.0)
    end
    achi2 = ROSESoss.chi2(mopt, ampobs, gs)
    rchi2 = (achi2)/(ROSE.nsamples(ampobs)-(logj.transform.dimension+(length(gs)-length(stations))))
    rachi2 = (achi2)/(ROSE.nsamples(ampobs)-(logj.transform.dimension+(length(gs)-length(stations))))
    machi2 = (achi2)/(ROSE.nsamples(ampobs))
    @info "Reduced chi square: $rchi2"

    #write summary files
    if plot_results
        samples = Array(echain)

        params = [NamedTuple{pnames}(samples[i,:]) for i in rand(1:size(samples,1), 500)]
        ims = ROSESoss.make_mean(logj.model, params)
        print(ims)
        p = ROSESoss.plot_im(ims)
        savefig(joinpath(@__DIR__, dir,  "MAP","mean_fit_$jscan.png"))

        #Plot some fit statistics stuff for MAP
        p1=plot_amp_comp(ampobs, gs, mopt)
        title!(p1,"Time: $tstart   χ²ᵣ = $rachi2")
        savefig(p1,joinpath(@__DIR__, dir,  "Residual", "amp_residuals_scan-$jscan.png"))

    end

    return merge((time = tstart, rchi2=rchi2, rampchi2=rachi2,
                    mampchi2=machi2,
                    logz=logz, logzerr=logzerr), opt)
end

function fit_scan(dir::String,
                              model,
                              visobs::ROSE.EHTObservation{F,V},
                              jscan::Int, sampler;
                              fitgains=false,
                              plot_results=true,
                              obschar=ROSESoss.ObsChar(),
                              kwargs...
                              ) where {F,V<:ROSE.EHTVisibilityDatum}

    tstart = round(visobs[1].time, digits=2)

    logj = ROSESoss.create_joint(model, visobs; fitgains=fitgains)

    chain, state, pnames = sampler(logj; kwargs...)
    logz = chain.logevidence
    logzerr = state.logzerr

    ess = 1/sum(abs2, chain[:weights])
    echain = sample(chain, Weights(vec(chain[:weights])), 2*Int(floor(ess)))
    @info "Nested sampling finished ess = $(ess)"

    DataFrame(echain) |> CSV.write(joinpath(@__DIR__,dir, "Chain",  "nested_chain_scan-$jscan.csv"))

    # Construct the image
    opt = NamedTuple{pnames}(Array(chain[end]))
    mopt = Soss.predict(logj.model, opt)[:img]
    intensitymap(mopt, 120,120,120,120)

    #Find a chi-square because yes...
    if fitgains
        ga = (AP=opt[:aAP], AZ=opt[:aAZ], JC=opt[:aJC], SM=opt[:aSM], AA=opt[:aAA],
                  LM=opt[:aLM], SP=opt[:aSP])
        gp = (AP=opt[:pAP], AZ=opt[:pAZ], JC=opt[:pJC], SM=opt[:pSM], AA=opt[:pAA],
                  LM=opt[:pLM], SP=opt[:pSP])
    else
        ga = (AP=1.0, AZ=1.0, JC=1.0, SM=1.0, AA=1.0,
                  LM=1.0, SP=1.0)
        gp = (AP=0.0, AZ=0.0, JC=0.0, SM=0.0, AA=0.0,
                  LM=0.0, SP=0.0)
    end
    vchi2 = ROSESoss.chi2(mopt, visobs, ga, gp)
    rchi2 = vchi2/(2*ROSE.nsamples(visobs)-logj.transform.dimension)
    mchi2 = (vchi2)/(2*ROSE.nsamples(visobs))
    @info "Reduced chi square: $rchi2"

    #write summary files
    if plot_results
        samples = Array(echain)

        params = [NamedTuple{pnames}(samples[i,:]) for i in rand(1:size(samples,1), 500)]
        ims = ROSESoss.make_mean(logj.model, params)
        p = ROSESoss.plot_im(ims)
        savefig(joinpath(@__DIR__, dir,  "MAP","mean_fit_$jscan.png"))

        #Plot some fit statistics stuff for MAP
        p1=plot_vis_comp(visobs, ga, gp, mopt)
        title!(p1,"Time: $tstart   χ²ᵣ = $rchi2")
        savefig(p1,joinpath(@__DIR__, dir,  "Residual", "amp_residuals_scan-$jscan.png"))

        Plots.closeall()
    end

    return merge((time = tstart, rchi2=rchi2,
                    mchi2=mchi2,
                    logz=logz, logzerr=logzerr), opt)
end



function hdf5_chain(dir, outname)
    dfchain, dfsum = load_chains(dir)
    h5open(outname, "w") do fid
        fid["time"] = dfsum[:,:time]
        fid["logz"] = dfsum[:,:logz]
        pid = create_group(fid, "params")
        for i in 1:length(dfchain)
            sid = create_group(pid, "scan$i")
            keys = names(dfchain[i])
            for k in keys
                write(sid, k, dfchain[i][:,Symbol(k)])
            end
        end
    end
end
