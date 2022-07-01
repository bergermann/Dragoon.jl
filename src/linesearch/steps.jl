#determines step length and/or normalization of step vector

function stepNorm(p,α,booster,hist,freqs,objFunction,mode; showtrace=false)
    #static, normalize p
    if mode == "unit"
        p[:] = p/pNorm(p)
    elseif mode == "max"
        p[:] = p/maximum(abs.(p))
    end
end
