
export tempLinear, TempLinear,
    tempExponential, TempExp




"""
    tempLinear(τ,(τ0,n)::Tuple{Real,Int})

Linear temperature decline by steps of `τ0/n`. Prevents negative temperatures.
"""
function tempLinear(τ,(τ0,n)::Tuple{Real,Int})
    τ_ = τ-τ0/n
    
    τ_ == 0 && (@warn "Warning: Temperature would go to zero!")
    τ_ < 0 && (@warn "Warning: Temperature would go below zero!")

    return max(τ_,0.)
end


"""
    TempLinear(τ0,n)

Callback for linear temperature manager option in [`simulatedAnnealing`](@ref). See 
[`tempLinear`](@ref).
"""
TempLinear(τ0,n) = Callback(tempLinear,(τ0,n))



function tempExponential(τ,(τ0,α,τmin)::Tuple{Real,Real})
    @assert 0 < α < 1 "α should be between 0 and 1."

    return max(α*τ,τmin)
end

"""
    TempExp(α,τmin=0)

Callback for exponential temperature manager option in [`simulatedAnnealing`](@ref). See 
[`tempExponential`](@ref).
"""
TempExp(τ0,α,τmin=0) = Callback(tempExponential,(τ0,α,τmin))


