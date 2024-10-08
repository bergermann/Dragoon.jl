###     how to calculate boost and objective functions

export getObjAna1d, getObjRef1d
export ObjAnalytical
export ObjRef1d, ObjRef1dLin, ObjRef1dSquare, ObjRef1dExp
export ObjRef1dTest
export ObjRef1dS, ObjRef1dSGD



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
    ref1 = normalize_range(abs.(getRef1d(booster,freqs)))
    ref2 = normalize_range(abs.(ref_goal))

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
    
    ref1 = normalize_range(abs.(getRef1d(booster,freqs)))
    ref2 = normalize_range(abs.(ref_goal))

    gd1 = groupdelay(ref1,freqs)
    gd2 = groupdelay(ref2,freqs)

    return sum(scaling.(abs.(ref1-ref2)))*sum(scaling_gd.(abs.(gd1-gd2)))
end

"""
    ObjRef1dSGD(ref,scaling,scaling_gd)

Callback for calculating objective value with scaled analytical 1d reflectivity.
See [`getObjRef1dSDG`](@ref).
"""
ObjRef1dSGD(ref,scaling,scaling_gd) = Callback(getObjRef1dSDG,(ref,scaling,scaling_gd))


"""
    ObjRef1dSGD(ref,scaling)

Callback for calculating objective value with scaled analytical 1d reflectivity.
See [`getObjRef1dSDG`](@ref).
"""
ObjRef1dSGD(ref,scaling) = Callback(getObjRef1dSDG,(ref,scaling,scaling))



"""
    getObjRef1dSP(booster::Booster,freqs::Vector{Float64},
        (ref_goal,scaling,scaling_gd)::Tuple{Vector{ComplexF64},Function,Function})

Return objective value by analytical1d reflectivity ``(∑ scaling(||ref|-ref_goal||))*(∑ scaling_p(||p|-p_goal||))``
for the given `freqs` and `booster` state, where ref and ref_goal get scaled to their negative peak heights.
Takes difference in phase into account.
See[`ref1d`](@ref).
"""
function getObjRef1dSP(booster::Booster,freqs::Vector{Float64},
        (ref_goal,scaling,scaling_p)::Tuple{Vector{ComplexF64},Function,Function})
    
    ref = getRef1d(booster,freqs)
    ref1 = normalize_range(abs.(ref))
    ref2 = normalize_range(abs.(ref_goal))

    ϕ1 = angle.(ref).-angle(ref)
    ϕ2 = angle.(ref_goal).-angle(ref_goal)


    return sum(scaling.(abs.(ref1-ref2)))*sum(scaling_p.(abs.(ϕ1-ϕ2)))
end

"""
    ObjRef1dSP(ref,scaling,scaling_p)

Callback for calculating objective value with scaled analytical 1d reflectivity.
See [`getObjRef1dSP`](@ref).
"""
ObjRef1dSP(ref,scaling,scaling_p) = Callback(getObjRef1dSP,(ref,scaling,scaling_p))


"""
    ObjRef1dSP(ref,scaling)

Callback for calculating objective value with scaled analytical 1d reflectivity.
See [`getObjRef1dSP`](@ref).
"""
ObjRef1dSP(ref,scaling) = Callback(getObjRef1dSP,(ref,scaling,scaling))
