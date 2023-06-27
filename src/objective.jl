###     how to calculate boost and objective functions

export getObjAna1d
export ObjAnalytical, ObjRef1dLin, ObjRef1dSquare, ObjRef1dExp

function getObjAna1d(booster::Booster,freqs::Array{Float64},args::Tuple{})
    return -minimum(boost1d(pos2dist(booster.pos; thickness=booster.thickness),
        freqs; eps=booster.epsilon,thickness=booster.thickness))
end

const ObjAnalytical = Callback(getObjAna1d)

# args = (ref_goal,)
function getObjRef1d(booster::Booster,freqs::Array{Float64},args::Tuple{Array{ComplexF64}})
    return sum(abs.((ref1d(pos2dist(booster.pos; thickness=booster.thickness),
        freqs; eps=booster.epsilon,thickness=booster.thickness)-args[1])))
end

# args = (ref_goal,scale_func)
function getObjRef1d(booster::Booster,freqs::Array{Float64},args::Tuple{Array{ComplexF64},Function})
    scale = args[2]

    return sum(scale.(abs.((ref1d(pos2dist(booster.pos; thickness=booster.thickness),
        freqs; eps=booster.epsilon,thickness=booster.thickness)-args[1]))))
end

const ObjRef1dLin(ref0) = Callback(getObjRef1d,(ref0,))
const ObjRef1dSquare(ref0) = Callback(getObjRef1d,(ref0,x->x^2))
const ObjRef1dExp(ref0) = Callback(getObjRef1d,(ref0,exp))