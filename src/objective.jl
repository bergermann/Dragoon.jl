export getObjAna1d

#how to calculate boost and objective functions

function getObjAna1d(booster::Booster,freqs::Array{Float64})
    return -minimum(boost1d(pos2dist(booster.pos; thickness=booster.thickness),
        freqs; eps=booster.epsilon,thickness=booster.thickness))
end
