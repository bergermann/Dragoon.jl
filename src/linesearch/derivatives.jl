#calculates all necessary derivates for solver

function firstDerivative(g,h,booster,hist,freqs,objFunction,Δx,mode)
    updateHist!(booster,hist,freqs,objFunction)
    moveCommand(booster,[(1,Δx)])
    if mode == "double"
        for i in 1:booster.ndisk
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(i,-2Δx)])
            updateHist!(booster,hist,freqs,objFunction)

            g[i] = (hist[2].objvalue-hist[1].objvalue)/2Δx

            if i != booster.ndisk
                moveCommand(booster,[(i,Δx),(i+1,Δx)])
            else
                moveCommand(booster,[(i,Δx)])
            end
        end
    else
        for i in 1:booster.ndisk
            updateHist!(booster,hist,freqs,objFunction)

            g[i] = (hist[1].objvalue-hist[i+1].objvalue)

            if i != booster.ndisk
                moveCommand(booster,[(i,-Δx),(i+1,Δx)])
            end
        end
    end
end

function secondDerivative(g,h,booster,hist,freqs,objFunction,Δx1,Δx2,mode)
    updateHist!(booster,hist,freqs,objFunction)
    moveCommand(booster,[(1,Δx1)])
    if mode == "double"
        for i in 1:booster.ndisk
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(i,-2Δx1)])
            updateHist!(booster,hist,freqs,objFunction)

            g[i] = (hist[2].objvalue-hist[1].objvalue)/2Δx1

            if i != booster.ndisk
                moveCommand(booster,[(i,Δx1),(i+1,Δx1)])
            else
                moveCommand(booster,[(i,Δx1)])
            end
        end
        for i in 1:booster.ndisk, j in 1:booster.ndisk
            moveCommand(booster,[(i,Δx2),(j,Δx2)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(j,-Δx2)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(i,-Δx2),(j,Δx2)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(j,-Δx2)])
            updateHist!(booster,hist,freqs,objFunction)
            h[i,j] = (hist[4].objvalue-hist[3].objvalue-
                        hist[2].objvalue+hist[1].objvalue)/Δx2^2
        end
    end
end
