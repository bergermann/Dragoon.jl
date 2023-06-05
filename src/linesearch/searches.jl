###     handles the step and stopping

export searchStandard, searchExtSteps, searchExtDist
export SearchStandard, SearchExtendedSteps, SearchExtendedDist, SearchTrueNewton

# args = (ϵls,kmax)
function searchStandard(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    k = 0

    while k < args[2]
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue - hist[2].objvalue > -args[1]
            move(booster,-α*p; additive=true)
            updateHist!(booster,hist,freqs,objFunction)

            break
        end
        k += 1
    end

    return k
end

const SearchStandard(ϵls,kmax) = Callback(searchStandard,(ϵls,kmax))



# args = (kmax,)
function searchExtSteps(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    for k in 1:args[1]
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)
    end

    k0 = argmin((x->x.objvalue).(hist[1:(args[1]+1)]))

    move(booster,hist[k0].pos; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    return args[1]-k0+1
end

const SearchExtendedSteps(kmax) = Callback(searchExtSteps,(kmax,))



# args = (d,)
function searchExtDist(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    kmax = args[1]/(α*pNorm(p))

    for i in 1:kmax
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)
    end

    k0 = argmin((x->x.objvalue).(hist[1:(kmax+1)]))

    move(booster,hist[k0].pos; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    return kmax-k0+1
end

const SearchExtendedDist(d) = Callback(searchExtDist,(d,))



#args = ()
function searchTrueNewton(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    move(booster,α*p; additive=true)
end

const SearchTrueNewton = Callback(searchTrueNewton)