export searchStandard, searchExtSteps, searchExtDist

#handles the step and stopping

function searchStandard(p,α,booster,hist,freqs,objFunction,ϵls,kmax;
                                                                showtrace=false)
    k = 0
    while k < kmax
        moveCommand(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)
        if hist[1].objvalue - hist[2].objvalue > -ϵls
            moveCommand(booster,-α*p; additive=true)
            updateHist!(booster,hist,freqs,objFunction)
            break
        end
        k += 1
    end

    return k
end

function searchExtSteps(p,α,booster,hist,freqs,objFunction,kmax; showtrace=false)
    for k in 1:kmax
        moveCommand(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)
    end

    k0 = argmin((x->x.objvalue).(hist[1:(kmax+1)]))

    moveCommand(booster,hist[k0].pos; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    return kmax-k0+1
end

function searchExtDist(p,α,booster,hist,freqs,objFunction,d; showtrace=false)
    kmax = d/(α*pNorm(p))
    for i in 1:kmax
        moveCommand(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)
    end

    k0 = argmin((x->x.objvalue).(hist[1:(kmax+1)]))

    moveCommand(booster,hist[k0].pos; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    return kmax-k0+1
end

function searchTrueNewton(p,α,booster,hist,freqs,objFunction; showtrace=false)
    moveCommand(booster,α*p; additive=true)
end
