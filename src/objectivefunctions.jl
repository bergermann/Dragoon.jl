###     how to calculate boost and objective functions

export getObjAna1d, getObjRef1d
export ObjAnalytical
export ObjRef, ObjRefLin, ObjRefSquare, ObjRefSquare, ObjRefExp



"""
    getObjAna1d(booster::Booster,freqs::Vector{Float64},args::Tuple{})

Return objective value by analytical1d boost for the given `freqs` and `booster`
state.
"""
function getObjAna1d(booster::Booster,freqs::Vector{Float64},args::Tuple{})
    return -minimum(boost1d(pos2dist(booster.pos; 
        disk_thickness=booster.thickness), freqs; 
        eps=booster.epsilon,thickness=booster.thickness,
        tand=booster.tand))
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
for the given `freqs` and `booster` state. 
"""
function getObjRef1d(booster::Booster,freqs::Vector{Float64},
        (ref_goal,)::Tuple{Vector{ComplexF64}})
    return sum(abs.(ref1d(pos2dist(booster.pos;
        disk_thickness=booster.thickness),freqs;
        eps=booster.epsilon,thickness=booster.thickness,
        tand=booster.tand)-ref_goal))
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

Return objective value by analytical1d reflectivity
``∑ scaling(|ref-ref_goal|)`` for the given `freqs` and `booster` state.
"""
function getObjRef1d(booster::Booster,freqs::Vector{Float64},
        (ref_goal,scaling)::Tuple{Vector{ComplexF64},Function})

    return sum(scaling.(abs.(ref1d(
        pos2dist(booster.pos; disk_thickness=booster.thickness), freqs;
        eps=booster.epsilon,thickness=booster.thickness,tand=booster.tand)
        -ref_goal)))
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