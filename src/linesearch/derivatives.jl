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
    secondDerivative(g,h,booster,hist,freqs,objFunction,args)

Calculate first and second derivative of objective function with finite differences and
write to `g`, `h`.

# args
- `Δx1`: Finite difference distance for first derivative.
- `Δx2`: Finite difference distance for second derivative.
- `mode`: Finite difference approach, one-side "single" or double-side "double".
"""
function secondDerivative(g,h,booster,hist,freqs,objFunction,(Δx1,Δx2,mode))
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

        for i in 1:booster.ndisk, j in 1:booster.ndisk
            if i == j
                move(booster,[(i,Δx2)])
                updateHist!(booster,hist,freqs,objFunction; force=true)
                
                move(booster,[(i,-Δx2)])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                move(booster,[(i,-Δx2)])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                h[i,i] = (hist[3].objvalue-2*hist[2].objvalue+hist[1].objvalue)/(Δx2^2)

                move(booster,[(i,Δx2)])
            else
                # x + h*e_i + h*e_j
                move(booster,[(i,Δx2),(j,Δx2)])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                # x + h*e_i - h*e_j
                move(booster,[(j,-2*Δx2)])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                # x - h*e_i + h*e_j
                move(booster,[(i,-2*Δx2),(j,2*Δx2)])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                # x - h*e_i - h*e_j
                move(booster,[(j,-2Δx2)])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                h[i,j] = (hist[4].objvalue-hist[3].objvalue-
                            hist[2].objvalue+hist[1].objvalue)/(4*Δx2^2)

                # h[i,j] = h[j,i] = (hist[4].objvalue-hist[3].objvalue-
                #             hist[2].objvalue+hist[1].objvalue)/(4*args[2]^2)

                move(booster,[(i,Δx2),(j,Δx2)])
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

"""
    Derivator2(Δx1,Δx2,mode)

Callback for second order derivator option in [`linesearch`](@ref). See 
[`secondDerivative`](@ref).
"""
const Derivator2(Δx1,Δx2,mode) = Callback(secondDerivative,(Δx1,Δx2,mode))





# args = (Δx1,Δx2,mode)
function secondDerivative_(g,h,booster,hist,freqs,objFunction,args)
    updateHist!(booster,hist,freqs,objFunction; force=true)

    p0 = copy(booster.pos)

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
                # move(booster,[(i,args[1])])
                move(booster,p0; additive=false)
            end
        end

        for i in 1:booster.ndisk, j in 1:booster.ndisk
            if i == j
                move(booster,[(i,args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)
                
                move(booster,[(i,-args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                move(booster,[(i,-args[2])])
                updateHist!(booster,hist,freqs,objFunction; force=true)

                h[i,i] = (hist[3].objvalue-2*hist[2].objvalue+hist[1].objvalue)/(args[2]^2)

                move(booster,[(i,args[2])])
                move(booster,p0; additive=false)
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

                # h[i,j] = (hist[4].objvalue-hist[3].objvalue-
                #             hist[2].objvalue+hist[1].objvalue)/(4*args[2]^2)

                h[i,j] = h[j,i] = (hist[4].objvalue-hist[3].objvalue-
                            hist[2].objvalue+hist[1].objvalue)/(4*args[2]^2)

                move(booster,[(i,args[2]),(j,args[2])])
                move(booster,p0; additive=false)
            end
        end
    else
        0
    end
end

const Derivator2_(Δx1,Δx2,mode) = Callback(secondDerivative_,(Δx1,Δx2,mode))



function secondDerivative__(g,h,booster,hist,freqs,objFunction,(Δx1,Δx2,mode))
    fx = updateHist!(booster,hist,freqs,objFunction; force=true)
    x0 = copy(booster.pos)
    
    if mode == "double"
        for i in 1:booster.ndisk
            move(booster,[(i,-Δx1)])
            updateHist!(booster,hist,freqs,objFunction)

            move(booster,[(i,Δx1)])
            updateHist!(booster,hist,freqs,objFunction)

            g[i] = (hist[1].objvalue-hist[2].objvalue)/2Δx1
            h[i,i] = (hist[1].objvalue+hist[2].objvalue-2fx)/Δx1^2

            for j in i+1:booster.ndisk
                move(booster,x0+e(booster.ndisk,(i,j)); additive=false)
                updateHist!(booster,hist,freqs,objFunction) #++

                move(booster,[(i,-2Δx2)])
                updateHist!(booster,hist,freqs,objFunction) #-+

                move(booster,[(j,-2Δx2)])
                updateHist!(booster,hist,freqs,objFunction) #--

                move(booster,[(i,2Δx2)])
                updateHist!(booster,hist,freqs,objFunction) #+-

                h[i,j] = h[j,i] = (hist[4].objvalue-hist[1].objvalue-hist[3].objvalue+
                    hist[2].objvalue)/4Δx^2
            end
        end
    else

    end

    move(booster,x0; additive=false)
end