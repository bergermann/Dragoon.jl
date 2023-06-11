###     how to calculate boost and objective functions

export getObjAna1d
export ObjAnalytical

function getObjAna1d(booster::Booster,freqs::Array{Float64},args::Tuple{})
    return -minimum(boost1d(pos2dist(booster.pos; thickness=booster.thickness),
        freqs; eps=booster.epsilon,thickness=booster.thickness))
end

const ObjAnalytical = Callback(getObjAna1d)
