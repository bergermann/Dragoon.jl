###     how to calculate boost and objective functions

export getObjAna1d, getObjRef1d
export ObjAnalytical
export ObjRef1d, ObjRef1dLin, ObjRef1dSquare, ObjRef1dExp
export ObjRef1dTest



"""
    getObjAna1d(booster::Booster,freqs::Vector{Float64},args::Tuple{})

Return objective value by analytical1d boost for the given `freqs` and `booster`
state. See[`boost1d`](@ref).
"""
function getObjAna1d(booster::Booster,freqs::Vector{Float64},args::Tuple{})
    return -minimum(getBoost1d(booster,freqs))
end

"""
    ObjAnalytical

Callback for calculating objective value with analytical1d boost.
See [`getObjAna1d`](@ref).
"""
ObjAnalytical = Callback(getObjAna1d)



"""
    getObjRef1d(booster::Booster,freqs::Vector{Float64},
        (ref_goal,)::Tuple{Vector{ComplexF64}})

Return objective value by analytical1d reflectivity ``∑|ref-ref_goal|``
for the given `freqs` and `booster` state. See[`ref1d`](@ref).
"""
function getObjRef1d(booster::Booster,freqs::Vector{Float64},
        (ref_goal,)::Tuple{Vector{ComplexF64}})
    return sum(abs.(getRef1d(booster,freqs)-ref_goal))
end

"""
    ObjRef1dLin(ref0)

Callback for calculating objective value with analytical1d reflectivity and
linear scaling. See [`getObjRef1d`](@ref).
"""
ObjRef1dLin(ref0) = Callback(getObjRef1d,(ref0,))



"""
    getObjRef1d(booster::Booster,freqs::Vector{Float64},
        (ref_goal,scaling)::Tuple{Vector{ComplexF64},Function})

Return objective value by analytical1d reflectivity ``∑ scaling(|ref-ref_goal|)``
for the given `freqs` and `booster` state. See[`ref1d`](@ref).
"""
function getObjRef1d(booster::Booster,freqs::Vector{Float64},
        (ref_goal,scaling)::Tuple{Vector{ComplexF64},Function})

    return sum(scaling.(abs.(getRef1d(booster,freqs)-ref_goal)))
end

"""
    ObjRef1d(ref0,f)

Callback for calculating objective value with analytical1d reflectivity and
linear scaling. See [`getObjRef1d`](@ref).
"""
ObjRef1d(ref0,f) = Callback(getObjRef1d,(ref0,f))

"""
    ObjRef1dSquare(ref0)

Callback for calculating objective value with analytical1d reflectivity and
quadratic scaling. See [`getObjRef1d`](@ref).
"""
ObjRef1dSquare(ref0) = Callback(getObjRef1d,(ref0,x->x^2))

"""
    ObjRef1dExp(ref0)

Callback for calculating objective value with analytical1d reflectivity and
exponential scaling. See [`getObjRef1d`](@ref).
"""
ObjRef1dExp(ref0) = Callback(getObjRef1d,(ref0,exp))




"""
    getObj3d(booster::Booster,freqs::Vector{Float64},
        (M,L,gridwidth,dx)::Tuple{Int,Int,Real,Real})


Return objective value by 3d boost using Bessel-modes for the given `freqs` and `booster`
state. See [`boost3d`](@ref).
"""
function getObj3d(booster::Booster,freqs::Vector{Float64},
        (M,L,gridwidth,dx)::Tuple{Int,Int,Real,Real})

    return -minimum(getBoost3d(booster,freqs,(M,L,gridwidth,dx)))
end

"""
    Obj3d(M,L,gridwidth,dx)

Callback for calculating objective value with 3d Bessel-mode boost.
See [`getObj3d`](@ref).
"""
Obj3d(M,L,gridwidth,dx) = Callback(getObj3d,(M,L,gridwidth,dx))



"""
    getObjRef1dTest(booster::Booster,freqs::Vector{Float64},args::Tuple{})

Return objective value by minimum analytical1d reflectivity for the given `freqs` and
`booster` state. See[`boost1d`](@ref).
"""
function  getObjRef1dTest(booster::Booster,freqs::Vector{Float64},args::Tuple{})
    return maximum(abs.(getRef1d(booster,freqs)))
end


"""
    ObjRef1dTest()

Callback for calculating objective value with analytical 1d reflectivity.
See [`getObjRef1dTest`](@ref).
"""
ObjRef1dTest() = Callback(getObjRef1dTest,())



"""
    getObjRef1dS(booster::Booster,freqs::Vector{Float64},
        (ref_goal,scaling)::Tuple{Vector{ComplexF64},Function})

Return objective value by analytical1d reflectivity ``∑ scaling(||ref|-|ref_goal||)``
for the given `freqs` and `booster` state, where ref and ref_goal get scaled to their negative peak heights.
See[`ref1d`](@ref).
"""
function getObjRef1dS(booster::Booster,freqs::Vector{Float64},(ref_goal,scaling)::Tuple{Vector{ComplexF64},Function})
    ref1 = abs.(getRef1d(booster,freqs))
    ref1 /= minimum(ref1)
    ref2 = abs.(ref_goal)
    ref2 /= minimum(ref)

    return sum(@. scaling(abs(ref1-ref2)))
end

"""
    ObjRef1dS(ref,scaling)

Callback for calculating objective value with scaled analytical 1d reflectivity.
See [`getObjRef1dS`](@ref).
"""
ObjRef1dS(ref,scaling) = Callback(getObjRef1dS,(ref,scaling))



"""
    getObjRef1dSGD(booster::Booster,freqs::Vector{Float64},
        (ref_goal,scaling,scaling_gd)::Tuple{Vector{ComplexF64},Function,Function})

Return objective value by analytical1d reflectivity ``(∑ scaling(||ref|-ref_goal||))*(∑ scaling_gd(||gd|-gd_goal||))``
for the given `freqs` and `booster` state, where ref and ref_goal get scaled to their negative peak heights.
Takes difference in groupdelay into account.
See[`ref1d`](@ref).
"""
function getObjRef1dSGD(booster::Booster,freqs::Vector{Float64},
        (ref_goal,scaling,scaling_gd)::Tuple{Vector{ComplexF64},Function,Function})
    
    ref = getRef1d(booster,freqs)
    ref1 = abs.(ref)
    ref1 /= minimum(ref1)
    ref2 = abs.(ref_goal)
    ref2 /= minimum(ref)

    gd1 = groupdelay(ref,freqs)
    gd2 = groupdelay(ref_goal,freqs)

    return sum(scaling.(abs.(ref1-ref2)))*sum(scaling_gd.(abs.(gd1-gd2)))
end

"""
    ObjRef1dSGD(ref,scaling,scaling_gd)

Callback for calculating objective value with scaled analytical 1d reflectivity.
See [`getObjRef1dS`](@ref).
"""
ObjRef1dSGD(ref,scaling,scaling_gd) = Callback(getObjRef1dS,(ref,scaling,scaling_gd))


"""
    ObjRef1dSGD(ref,scaling)

Callback for calculating objective value with scaled analytical 1d reflectivity.
See [`getObjRef1dS`](@ref).
"""
ObjRef1dSGD(ref,scaling) = Callback(getObjRef1dS,(ref,scaling,scaling))
