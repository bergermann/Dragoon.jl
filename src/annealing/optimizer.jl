

export simulatedAnnealing








function simulatedAnnealing(booster::Booster,hist::Vector{State},freqs::Array{Float64},
        τ::Array{Float64},rmax::Float64,
        objFunction::Callback,
        unstuckinator::Callback;
        maxiter::Integer=Int(1e2),
        nreset::Int64=0,
        nresetterm::Int64=0,
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        traceevery::Int=1,
        resettimer::Bool=true,
        returntimes::Bool=false)

    if hasproperty(booster,:startingtime) && resettimer
        showtrace && println("Resetting starting time.")
        
        booster.startingtime = unow()
        booster.timestamp = unow()
    elseif hasproperty(booster,:codetimestamp)
        t0 = unow()

        if resettimer
            booster.timestamp = DateTime(0)
        end
    end

    trace = Vector{SATrace}(undef,round(Int,maxiter/traceevery)+1)

    updateHist!(booster,hist,freqs,objFunction)
    
    x = copy(booster.pos)
    objx = hist[1].objvalue
    
    xsol = copy(booster.pos)
    objsol = hist[1].objvalue

    iter = 0
    i = 0
    n_τ = length(τ)
    resetcounter = 0
    resetcounterterm = 0

    while iter < maxiter && i < n_τ
        if Int(iter%traceevery)==0
            trace[Int(iter/traceevery)+1] = SATrace(x,objx,xsol,objsol,τ[i+1],iter,
                                    booster.timestamp,booster.summeddistance)
        end

        iter += 1
        i += 1

        # updateHist!(booster,hist,freqs,objFunction; force=true)
        move(booster,x+findNeighbour(booster,rmax); additive=false)
        updateHist!(booster,hist,freqs,objFunction; force=true)

        if hist[1].objvalue <= objx || rand() <= thermal(objx,hist[1].objvalue,τ[i])
            x = copy(booster.pos)
            objx = hist[1].objvalue
        end

        if hist[1].objvalue <= objsol
            xsol = copy(booster.pos)
            objsol = hist[1].objvalue

            resetcounter = 0
            resetcounterterm = 0
        else
            resetcounter += 1
        end

        if nreset > 0 && resetcounter >= nreset
            showtrace && println("Resetting to best solution.")
            resetcounter = 0
            resetcounterterm += 1

            if nresetterm > 0 && resetcounterterm >= nresetterm
                showtrace && println("$nresetterm times resetted. Terminating.")

                break
            end
        end
    

        showtrace && i%showevery == 0 && printSAIter(booster,objx,objsol,τ[i],iter)

        ## unstucking, alter i
    end



    move(booster,xsol; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    if hasproperty(booster,:codetimestamp)
        if resettimer
            booster.codetimestamp = DateTime(0)
        end

        booster.codetimestamp += unow()-t0
    end

    term = printTermination(booster,hist,i,maxiter,showtrace)

    if returntimes
        return trace[1:min(round(Int,iter/traceevery)+1,length(trace))], term
    end

    return trace[1:min(round(Int,iter/traceevery)+1,length(trace))]
    # return trace
end
