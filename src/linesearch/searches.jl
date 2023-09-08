###     handles the step and stopping

export searchStandard, searchExtSteps, searchExtDist
export SearchStandard, SearchExtendedSteps, SearchExtendedDist, SearchTrueNewton

# args = (ϵls,kmax)
function searchStandard(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    updateHist!(booster,hist,freqs,objFunction)
    
    k = 0

    while k < args[2]
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue - hist[2].objvalue > -args[1]
            move(booster,α*p; additive=true)
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
    updateHist!(booster,hist,freqs,objFunction)
    
    k0 = 0
    obj0 = hist[1].objvalue
    pos0 = copy(booster.pos)

    for k in 1:args[1]
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue < obj0
            k0 = k
            obj0 = hist[1].objvalue
            pos0 = copy(booster.pos)
        end
    end

    move(booster,pos0; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    return k0
end

const SearchExtendedSteps(kmax) = Callback(searchExtSteps,(kmax,))



# args = (d,)
function searchExtDist(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    kmax = args[1]/(α*pNorm(p))

    updateHist!(booster,hist,freqs,objFunction)
    
    k0 = 0
    obj0 = hist[1].objvalue
    pos0 = copy(booster.pos)

    for i in 1:kmax
        move(booster,α*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue < obj0
            k0 = k
            obj0 = hist[1].objvalue
            pos0 = copy(booster.pos)
        end
    end

    # k0 = argmin((x->x.objvalue).(hist[1:(kmax+1)]))

    move(booster,pos0; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    return k0
end

const SearchExtendedDist(d) = Callback(searchExtDist,(d,))



#args = ()
function searchTrueNewton(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    move(booster,α*p; additive=true)
end

const SearchTrueNewton = Callback(searchTrueNewton)