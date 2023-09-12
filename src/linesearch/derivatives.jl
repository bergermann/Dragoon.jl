###     calculates all necessary derivates for solver

export firstDerivative, secondDerivative
export Derivator1, Derivator2

# args = (Δx,mode)
function firstDerivative(g,h,booster,hist,freqs,objFunction,args)
    updateHist!(booster,hist,freqs,objFunction; force=true)

    move(booster,[(1,args[1])])

    if args[2] == "double"
        for i in 1:booster.ndisk
            updateHist!(booster,hist,freqs,objFunction; force=true)
            move(booster,[(i,-2*args[1])])
            updateHist!(booster,hist,freqs,objFunction; force=true)

            g[i] = (hist[2].objvalue-hist[1].objvalue)/(2*args[1])

            if i != booster.ndisk
                move(booster,[(i,args[1]),(i+1,args[1])])
            else
                move(booster,[(i,args[1])])
            end
        end
    else
        for i in 1:booster.ndisk
            updateHist!(booster,hist,freqs,objFunction; force=true)

            g[i] = (hist[1].objvalue-hist[i+1].objvalue)/(args[1])

            if i != booster.ndisk
                move(booster,[(i,-args[1]),(i+1,args[1])])
            end
        end
    end
end

const Derivator1(Δx,mode) = Callback(firstDerivative,(Δx,mode))



# args = (Δx1,Δx2,mode)
function secondDerivative(g,h,booster,hist,freqs,objFunction,args)
    updateHist!(booster,hist,freqs,objFunction; force=true)

    move(booster,[(1,args[1])])
    
    if args[3] == "double"
        for i in 1:booster.ndisk
            updateHist!(booster,hist,freqs,objFunction; force=true)
            move(booster,[(i,-2*args[1])])
            updateHist!(booster,hist,freqs,objFunction; force=true)

            g[i] = (hist[2].objvalue-hist[1].objvalue)/(2*args[1])

            if i != booster.ndisk
                move(booster,[(i,args[1]),(i+1,args[1])])
            else
                move(booster,[(i,args[1])])
            end
        end

        for i in 1:booster.ndisk, j in 1:i
            if i == j
                move(booster,[(i,args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)
                
                move(booster,[(i,-args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                move(booster,[(i,-args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                h[i,i] = (hist[3].objvalue-2*hist[2].objvalue+hist[1].objvalue)/(args[2]^2)

                move(booster,[(i,args[2])])
            else
                # x + h*e_i + h*e_j
                move(booster,[(i,args[2]),(j,args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                # x + h*e_i - h*e_j
                move(booster,[(j,-2*args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                # x - h*e_i + h*e_j
                move(booster,[(i,-2*args[2]),(j,2*args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                # x - h*e_i - h*e_j
                move(booster,[(j,-2args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                h[i,j] = h[j,i] = (hist[4].objvalue-hist[3].objvalue-
                            hist[2].objvalue+hist[1].objvalue)/(4*args[2]^2)

                move(booster,[(i,args[2]),(j,args[2])])
            end
        end
    else
                # # x + h*e_i + h*e_j
        # move(booster,[(i,args[2]),(j,args[2])])
        # updateHist!(booster,hist,freqs,objFunction; force=true)

        # # x + h*e_i
        # move(booster,[(j,-args[2])])
        # updateHist!(booster,hist,freqs,objFunction; force=true)

        # # x + h*e_j
        # move(booster,[(i,-args[2]),(j,args[2])])
        # updateHist!(booster,hist,freqs,objFunction; force=true)

        # # x 
        # move(booster,[(j,-args[2])])
        # updateHist!(booster,hist,freqs,objFunction; force=true)
        
        # h[i,j] = (hist[4].objvalue-hist[3].objvalue-
        #             hist[2].objvalue+hist[1].objvalue)/(args[2]^2)
    end
end

const Derivator2(Δx1,Δx2,mode) = Callback(secondDerivative,(Δx1,Δx2,mode))