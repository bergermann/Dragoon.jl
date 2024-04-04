
export unstuckDont



"""
    unstuckDont(booster,hist,freqs,objFunction,i,args; showtrace=false)

Don't attempt to unstuck booster, terminate instead.
"""
function unstuckDont(booster,hist,freqs,objFunction,i,args; showtrace=false)
    showtrace && println("No unstucking tried. Terminating.")

    return true
end



"""
    unstuckTemp(booster,hist,freqs,objFunction,objsol,τ,(β,threshold); showtrace=false)

Attempt to unstuck booster by increasing temperature `τ` by factor `β`.
Don't unstuck if objective value is below `threshold`.
"""
function unstuckTemp(booster,hist,freqs,objFunction,objsol,τ,(β,threshold); showtrace=false)
    if objsol < threshold
        showtrace && println("Unstuck threshold reached. Terminating.")
        
        return true
    end

    @assert β > 1 "Unstuck factor β should be larger than 1."
    
    τ *= β
    
    showtrace && println("Unstuck successfull.")

    return false
end

"""
    UnstuckTemp(β,threshold)

Callback for unstucker that rescales the temperature in [`simulatedAnnealing`](@ref). See 
[`unstuckTemp`](@ref).
"""
UnstuckTemp(β,threshold)