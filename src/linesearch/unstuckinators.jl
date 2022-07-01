export unstuckDont, unstuckRandom, unstuckCoord

#tries to unstuck system if termination is unsatisfactory

function unstuckDont(booster,hist,freqs,objFunction; showtrace=false)
    showtrace && println("No unstucking tried. Terminating.")

    return true
end

#shift discs independently and uniform randomly within [-d,d] if objective
#threshold is not reached
function unstuckRandom(booster,hist,freqs,objFunction,d,threshold;
                                                                showtrace=false)
    if hist[1].objvalue > threshold
        moveCommand(booster,d*(2*rand(booster.ndisk).-1); additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        showtrace && println("Unstuck successfull.")

        return false
    else
        showtrace && println("Unstuck threshold reached. Terminating.")

        return true
    end
end

#search every coordinate axis for better points
#search along distance d divided into div steps
#enforce movement by setting dmin to > 0
function unstuckCoord(booster,hist,freqs,objFunction,d,dmin,div::Int,threshold;
                                showtrace=false)
    if hist[1].objvalue > threshold
        pos0 = copy(booster.pos)

        for i in 1:booster.ndisk
            if dmin > 0.
                moveCommand(booster,[(i,dmin)]; additive=true)
                updateHist!(booster,hist,freqs,objFunction)
            end

            δ = (d-dmin)/div

            for j in 1:div
                moveCommand(booster,[(i,δ)])
                updateHist!(booster,hist,freqs,objFunction)
            end

            moveCommand(booster,[(i,-(d+dmin))])
            updateHist!(booster,hist,freqs,objFunction)

            for j in 1:div
                moveCommand(booster,[(i,-δ)])
                updateHist!(booster,hist,freqs,objFunction)
            end
        end

        k0 = argmin((x->x.objvalue).(hist[1:(2*booster.ndisk*(div+1))]))

        if hist[k0].pos != pos0
            moveCommand(booster,hist[k0].pos; additive=false)
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
