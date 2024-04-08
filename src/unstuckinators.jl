###     tries to unstuck system if termination is unsatisfactory

export unstuckDont, unstuckRandom, unstuckCoord
export UnstuckDont, UnstuckRandom, UnstuckCoord



"""
    unstuckDont(booster,hist,freqs,objFunction,args; showtrace=false)

Don't attempt to unstuck booster, terminate instead.
"""
function unstuckDont(booster,hist,freqs,objFunction,args; showtrace=false)
    showtrace && println("No unstucking tried. Terminating.")

    return true
end

"""
    UnstuckDont

Callback for default unstuck option that doesn't unstuck.
See [`unstuckDont`](@ref).
"""
UnstuckDont = Callback(unstuckDont)



"""
    unstuckRandom(booster,hist,freqs,objFunction,(d,threshold); showtrace=false)

Attempt to unstuck booster by shifting every disc by a uniformly random value
``∈`` `[-d,d]`. Don't unstuck if objective value is below `threshold`.
"""
function unstuckRandom(booster,hist,freqs,objFunction,(d,threshold); showtrace=false)
    if hist[1].objvalue > threshold
        move(booster,d*(2*rand(booster.ndisk).-1); additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        showtrace && println("Unstuck successfull.")

        return false
    else
        showtrace && println("Unstuck threshold reached. Terminating.")

        return true
    end
end

"""
    UnstuckRandom(d,threshold)

Callback for random unstuck option. See [`unstuckRandom`](@ref).
"""
UnstuckRandom(d,threshold) = Callback(unstuckRandom,(d,threshold))



"""
    unstuckCoord(booster,hist,freqs,objFunction,
        (d,dmin,div,threshold); showtrace=false)

Attempt to unstuck booster by performing search for `div` steps from `dmin`
to `d` along every coordinate axis. Don't unstuck if objective value is below 
`threshold`.
"""
function unstuckCoord(booster,hist,freqs,objFunction,
        (d,dmin,div,threshold); showtrace=false)

    if hist[1].objvalue > threshold
        pos0 = copy(booster.pos)

        for i in 1:booster.ndisk
            if dmin > 0.
                move(booster,[(i,dmin)]; additive=true)
                updateHist!(booster,hist,freqs,objFunction)
            end

            δ = (d-dmin)/div

            for _ in 1:div
                move(booster,[(i,δ)])
                updateHist!(booster,hist,freqs,objFunction)
            end

            move(booster,[(i,-(d+dmin))])
            updateHist!(booster,hist,freqs,objFunction)

            for _ in 1:div
                move(booster,[(i,-δ)])
                updateHist!(booster,hist,freqs,objFunction)
            end
        end

        k0 = argmin((x->x.objvalue).(hist[1:(2*booster.ndisk*(div+1))]))

        if hist[k0].pos != pos0
            move(booster,hist[k0].pos; additive=false)
            showtrace && println("Unstuck successfull.")

            return false
        else
            showtrace && println("Unstuck unsuccessfull. Terminating.")

            return true
        end
    else
        showtrace && println("Unstuck threshold reached. Terminating.")

        return true
    end
end


"""
    UnstuckCoord(d,dmin,div,threshold)

Callback for unstuck search along coordinate axes. See [`unstuckCoord`](@ref).
"""
UnstuckCoord(d,dmin,div,threshold) = 
    Callback(unstuckCoord,(d,dmin,div,threshold))
