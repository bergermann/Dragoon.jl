###     determines step length and/or normalization of step vector

export stepNorm


# args = (mode,)
function stepNorm(p,α,booster,hist,freqs,objFunction,args; showtrace=false)
    #static, normalize p
    if args[1] == "unit"
        p[:] = p/pNorm(p)
    elseif args[1] == "max"
        p[:] = p/maximum(abs.(p))
    end
end
