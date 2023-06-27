###     how to calculate boost and objective functions

export getObjAna1d
export ObjAnalytical

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

const ObjRefLin(ref0) = Callback(getObjRef1d,(ref0,))
const ObjRefSquare(ref0) = Callback(getObjRef1d,(ref0,x->x^2))
const ObjRefExp(ref0) = Callback(getObjRef1d,(ref0,exp))