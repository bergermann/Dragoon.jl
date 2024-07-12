###     how to calculate boost and objective functions

export getObjAna1d, getObjRef1d
export ObjAnalytical
export ObjRef, ObjRefLin, ObjRefSquare, ObjRefSquare, ObjRefExp
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
    ObjRefLin(ref0)

Callback for calculating objective value with analytical1d reflectivity and
linear scaling. See [`getObjRef1d`](@ref).
"""
ObjRefLin(ref0) = Callback(getObjRef1d,(ref0,))



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
    ObjRef(ref0,f)

Callback for calculating objective value with analytical1d reflectivity and
linear scaling. See [`getObjRef1d`](@ref).
"""
ObjRef(ref0,f) = Callback(getObjRef1d,(ref0,f))

"""
    ObjRefSquare(ref0)

Callback for calculating objective value with analytical1d reflectivity and
quadratic scaling. See [`getObjRef1d`](@ref).
"""
ObjRefSquare(ref0) = Callback(getObjRef1d,(ref0,x->x^2))

"""
    ObjRefExp(ref0)

Callback for calculating objective value with analytical1d reflectivity and
exponential scaling. See [`getObjRef1d`](@ref).
"""
ObjRefExp(ref0) = Callback(getObjRef1d,(ref0,exp))




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