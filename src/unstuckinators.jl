###     tries to unstuck system if termination is unsatisfactory

export unstuckDont, unstuckRandom, unstuckCoord
export UnstuckDont, UnstuckRandom, UnstuckCoord

function unstuckDont(booster,hist,freqs,objFunction,args; showtrace=false)
    showtrace && println("No unstucking tried. Terminating.")

    return true
end

const UnstuckDont = Callback(unstuckDont)

# shift discs independently and uniform randomly within [-d,d] if objective
# threshold is not reached
# args = (d,threshold)
function unstuckRandom(booster,hist,freqs,objFunction,args; showtrace=false)
    if hist[1].objvalue > args[2]
        move(booster,args[1]*(2*rand(booster.ndisk).-1); additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        showtrace && println("Unstuck successfull.")

        return false
    else
        showtrace && println("Unstuck threshold reached. Terminating.")

        return true
    end
end

UnstuckRandom(d,threshold) = Callback(unstuckRandom,(d,threshold))

# search every coordinate axis for better points
# search along distance d divided into div steps
# enforce movement by setting dmin to > 0
# args = (d,dmin,div::Int,threshold)
function unstuckCoord(booster,hist,freqs,objFunction,args; showtrace=false)
    threshold = -abs(args[2])

    if hist[1].objvalue > threshold
        pos0 = copy(booster.pos)

        for i in 1:booster.ndisk
            if args[2] > 0.
                move(booster,[(i,args[2])]; additive=true)
                updateHist!(booster,hist,freqs,objFunction)
            end

            δ = (args[1]-args[2])/args[3]

            for j in 1:args[3]
                move(booster,[(i,δ)])
                updateHist!(booster,hist,freqs,objFunction)
            end

            move(booster,[(i,-(args[1]+args[2]))])
            updateHist!(booster,hist,freqs,objFunction)

            for j in 1:args[3]
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

const UnstuckCoord(d,dmin,div,threshold) = Callback(unstuckCoord,(d,dmin,div,threshold))
