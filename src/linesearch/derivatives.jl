###     calculates all necessary derivates for solver

export firstDerivative, secondDerivative
export Derivator1, Derivator2

"""
    firstDerivative(g,h,booster,hist,freqs,objFunction,args)

Calculate first derivative of objective function with finite differences and write to `g`.

# args
- `Δx`: Finite difference distance.
- `mode`: Finite difference approach, one-side "single" or double-side "double".
"""
function firstDerivative(g,h,booster,hist,freqs,objFunction,(Δx,mode))
    updateHist!(booster,hist,freqs,objFunction; force=true)

    move(booster,[(1,Δx)])

    if mode == "double"
        for i in 1:booster.ndisk
            updateHist!(booster,hist,freqs,objFunction; force=true)
            move(booster,[(i,-2*Δx)])
            updateHist!(booster,hist,freqs,objFunction; force=true)

            g[i] = (hist[2].objvalue-hist[1].objvalue)/(2*Δx)

            if i != booster.ndisk
                move(booster,[(i,Δx),(i+1,Δx)])
            else
                move(booster,[(i,Δx)])
            end
        end
    else
        for i in 1:booster.ndisk
            updateHist!(booster,hist,freqs,objFunction; force=true)

            g[i] = (hist[1].objvalue-hist[i+1].objvalue)/(Δx)

            if i != booster.ndisk
                move(booster,[(i,-Δx),(i+1,Δx)])
            end
        end
    end
end

"""
    Derivator1(Δx,mode)

Callback for first order derivator option in [`linesearch`](@ref). See 
[`firstDerivative`](@ref).
"""
const Derivator1(Δx,mode) = Callback(firstDerivative,(Δx,mode))



"""
    secondDerivative(g,h,booster,hist,freqs,objFunction,(Δx1,Δx2,mode))

Calculate first and second derivative of objective function with finite differences and
write to `g`, `h`.

# args
- `Δx1`: Finite difference distance for first derivative.
- `Δx2`: Finite difference distance for second derivative.
- `mode`: Finite difference approach, one-side "single" or double-side "double".
"""
function secondDerivative(g,h,booster,hist,freqs,objFunction,(Δx1,Δx2,mode))
    fx = updateHist!(booster,hist,freqs,objFunction; force=true)
    x0 = copy(booster.pos)
    
    if mode == "double"
        for i in 1:booster.ndisk
            move(booster,x0-Δx1*e(booster.ndisk,[i,]); additive=false)
            updateHist!(booster,hist,freqs,objFunction)

            move(booster,[(i,2Δx1)])
            updateHist!(booster,hist,freqs,objFunction)

            g[i] = (hist[1].objvalue-hist[2].objvalue)/2Δx1
            h[i,i] = (hist[1].objvalue+hist[2].objvalue-2fx)/Δx1^2

            for j in i+1:booster.ndisk
                move(booster,x0+Δx2*e(booster.ndisk,[i,j]); additive=false)
                updateHist!(booster,hist,freqs,objFunction) #++

                move(booster,[(i,-2Δx2)])
                updateHist!(booster,hist,freqs,objFunction) #-+

                move(booster,[(j,-2Δx2)])
                updateHist!(booster,hist,freqs,objFunction) #--

                move(booster,[(i,2Δx2)])
                updateHist!(booster,hist,freqs,objFunction) #+-

                h[i,j] = h[j,i] = (hist[4].objvalue-hist[1].objvalue-hist[3].objvalue+
                    hist[2].objvalue)/4Δx2^2
            end
        end
    else
        for i in 1:booster.ndisk
            if Δx1 < Δx2
                move(booster,x0+Δx1*e(booster.ndisk,[i,]); additive=false)
                updateHist!(booster,hist,freqs,objFunction)

                g[i] = (hist[1].objvalue-fx)/Δx1

                move(booster,[i,Δx2-Δx1])
                fx_i = updateHist!(booster,hist,freqs,objFunction)

                move(booster,[i,Δx2])
                updateHist!(booster,hist,freqs,objFunction)

                h[i,i] = (hist[1].objvalue-2fx_i+fx)/Δx2^2
            elseif Δx2 <= Δx1 < 2Δx2
                move(booster,x0+Δx2*e(booster.ndisk,[i,]); additive=false)
                fx_i = updateHist!(booster,hist,freqs,objFunction)

                move(booster,[(i,Δx1-Δx2)])
                updateHist!(booster,hist,freqs,objFunction)

                g[i] = (hist[1].objvalue-fx)/Δx1

                move(booster,[(i,2Δx2-Δx1)])
                updateHist!(booster,hist,freqs,objFunction)

                h[i,i] = (hist[1].objvalue-2fx_i+fx)/Δx2^2
            else
                move(booster,x0+Δx2*e(booster.ndisk,[i,]); additive=false)
                fx_i = updateHist!(booster,hist,freqs,objFunction)

                move(booster,[(i,Δx2)])
                updateHist!(booster,hist,freqs,objFunction)

                h[i,i] = (hist[1].objvalue-2fx_i+fx)/Δx2^2

                move(booster,[(i,Δx1-2Δx2)])
                updateHist!(booster,hist,freqs,objFunction)

                g[i] = (hist[1].objvalue-fx)/Δx1
            end

            for j in i+1:booster.ndisk
                move(booster,x0+Δx2*e(booster.ndisk,[i,j]); additive=false)
                updateHist!(booster,hist,freqs,objFunction)

                move(booster,[(i,-Δx2)])
                updateHist!(booster,hist,freqs,objFunction)

                h[i,j] = h[j,i] = (hist[2].objvalue-hist[1].objvalue-fx_i-fx)/Δx2^2
            end
        end
    end

    move(booster,x0; additive=false)

    return
end

"""
    Derivator2(Δx1,Δx2,mode)

Callback for second order derivator option in [`linesearch`](@ref). See 
[`secondDerivative`](@ref).
"""
Derivator(Δx1,Δx2,mode) = Callback(secondDerivative,(Δx1,Δx2,mode))

