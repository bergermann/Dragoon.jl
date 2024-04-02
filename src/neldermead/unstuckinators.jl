



export unstuckDont
export unstuckNewSimplex, UnstuckNew
export unstuckExpandSimplex, UnstuckExpand


"""
    unstuckDont(booster,hist,freqs,objFunction,simplexObj,x,f,(); showtrace=false)

Default Nelder-Mead unstucker that just terminates.
"""
function unstuckDont(booster,hist,freqs,objFunction,simplexObj,x,f,(); showtrace=false)
    showtrace && println("No unstucking tried. Terminating.")

    return true
end



"""
    function unstuckNewSimplex(booster,hist,freqs,objFunction,simplexObj,x,f,
        (initSimplex,best,threshold); showtrace=false)

Create newly initialized simplex to unstuck Nelder-Mead.

# args
- initSimplex: Simplex initializer.
- best: Construct new simplex around best or worst old vertex.
- threshold: Objective threshold under which no further unstucking is tried.
"""
function unstuckNewSimplex(booster,hist,freqs,objFunction,simplexObj,x,f,
        (initSimplex,best,threshold); showtrace=false)

    if minimum(f) < threshold
        showtrace && println("Unstuck threshold reached. Terminating.")

        return true
    end

    best ? (x0 = x[:,argmin(f)]) : (x0 = x[:,argmax(f)])

    x[:] = initSimplex.func(x0,initSimplex.args)
    f[:] = simplexObj.func(x,collect(1:booster.ndisk+1),booster,hist,freqs,
                                            objFunction,simplexObj.args)

    showtrace && println("Unstuck successfull.")

    return false
end


"""
    UnstuckNew(init,best,threshold)

Callback for Nelder-Mead unstuck option by creating new simplex in [`nelderMead`](@ref).
See [`unstuckNewSimplex`](@ref).
"""
UnstuckNew(init,best,threshold) = Callback(unstuckNewSimplex,(init,best,threshold))




"""
    unstuckExpandSimplex(booster,hist,freqs,objFunction,simplexObj,x,f,
        (δ,threshold); showtrace=false)

Create soft new simplex by expand operation away from worst vertex.

# args
- δ: Expansion factor.
- threshold: Objective threshold under which no further unstucking is tried.
"""
function unstuckExpandSimplex(booster,hist,freqs,objFunction,simplexObj,x,f,
        (δ,threshold); showtrace=false)

    if minimum(f) < threshold
        showtrace && println("Unstuck threshold reached. Terminating.")
        
        return true
    end
    
    if δ < 1
        error("Unstucking requires expansion with δ > 1.")
    end
    
    for i in 1:booster.ndisk
        v = x[:,end]+δ*(x[:,i]-x[:,end])
        x[:,j] = v
    end

    f[1:end-1] = simplexObj.func(x,collect(1:booster.ndisk),booster,hist,freqs,objFunction,simplexObj.args)

    showtrace && println("Unstuck successfull.")

    return false
end

"""
    UnstuckExpand(δ,threshold)

Callback for Nelder-Mead unstuck option by expanding old simplex in [`nelderMead`](@ref).
See [`unstuckExpandSimplex`](@ref).
"""
const UnstuckExpand(δ,threshold) = Callback(unstuckExpandSimplex,(δ,threshold))
