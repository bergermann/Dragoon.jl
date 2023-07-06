###     how to calculate boost and objective functions

export getObjAna1d, getObjRef1d
export ObjAnalytical
export ObjRef, ObjRefLin, ObjRefSquare, ObjRefSquare, ObjRefExp

function getObjAna1d(booster::Booster,freqs::Vector{Float64},args::Tuple{})
    return -minimum(boost1d(pos2dist(booster.pos; thickness=booster.thickness),
        freqs; eps=booster.epsilon,thickness=booster.thickness))
end

const ObjAnalytical = Callback(getObjAna1d)

# args = (ref_goal,)
function getObjRef1d(booster::Booster,freqs::Vector{Float64},args::Tuple{Vector{Float64}})
    return sum(abs.(ref1d(pos2dist(booster.pos; thickness=booster.thickness),
        freqs; eps=booster.epsilon,thickness=booster.thickness)-args[1]))
end

const ObjRefLin(ref0) = Callback(getObjRef1d,(ref0,))

# args = (ref_goal,scaling)
function getObjRef1d(booster::Booster,freqs::Vector{Float64},args::Tuple{Vector{Float64},Function})
    scaling = args[2]

    return sum(scaling.(abs.(ref1d(pos2dist(booster.pos; thickness=booster.thickness),
        freqs; eps=booster.epsilon,thickness=booster.thickness)-args[1])))
end

const ObjRef(ref0,Function) = Callback(getObjRef1d,(ref0,Function))
const ObjRefSquare(ref0) = Callback(getObjRef1d,(ref0,x->x^2))
const ObjRefExp(ref0) = Callback(getObjRef1d,(ref0,exp))