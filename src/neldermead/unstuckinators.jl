



export unstuckDont
export unstuckNewSimplex, UnstuckNew
export unstuckExpandSimplex, UnstuckExpand



function unstuckDont(booster,hist,freqs,objFunction,simplexObj,x,f,args; showtrace=false)
    showtrace && println("No unstucking tried. Terminating.")

    return true
end

function unstuckNewSimplex(booster,hist,freqs,objFunction,simplexObj,x,f,(initSimplex,best,threshold); showtrace=false)
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

const UnstuckNew(init,best,threshold) = Callback(unstuckNewSimplex,(init,best,threshold))


function unstuckExpandSimplex(booster,hist,freqs,objFunction,simplexObj,x,f,(δ,threshold); showtrace=false)
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

const UnstuckExpand(δ,threshold) = Callback(unstuckExpandSimplex,(δ,threshold))
