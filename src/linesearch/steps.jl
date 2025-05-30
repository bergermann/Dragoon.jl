###     determines step length and/or normalization of step vector

export stepNorm
export StepNorm

"""
    stepNorm(p,α,booster,hist,freqs,objFunction,(mode,); showtrace=false)

Normalize descend direction p to unit or maximum entry.

# args
- `mode`: "unit" or "max".
"""
function stepNorm(p,α,booster,hist,freqs,objFunction,(mode,); showtrace=false)
    #static, normalize p
    if mode == "unit"
        p[:] = p/norm(p)
    elseif mode == "max"
        p[:] = p/maximum(abs.(p))
    end
end

"""
    StepNorm(mode)

Callback for step option in [`linesearch`](@ref). mode is either "unit" or "max". See
[`stepNorm`](@ref).
"""
StepNorm(mode) = Callback(stepNorm,(mode,))
