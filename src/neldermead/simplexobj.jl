###     functions to get objective values for every vertex

export getSimplexObj

function getSimplexObj(x::Matrix{Float64},
                        booster::Booster,
                        hist::Array{State},
                        freqs::Array{Float64},
                        objFunction::Callback,
                        args::Tuple{};
                        reset=false)
    
    reset && (xc = copy(booster.pos))

    for i in axes(x,2)
        move(booster,x[:,i]; additive=false)
        updateHist!(booster,hist,freqs,objFunction; force=true)
    end

    reset && move(booster,xc; additive=false)

    return (a->a.objvalue).(hist[axes(x,2)])
end

function getSimplexObj(x::Matrix{Float64},
                        indices::Vector{Int},
                        booster::Booster,
                        hist::Array{State},
                        freqs::Array{Float64},
                        objFunction::Callback,
                        args::Tuple{};
                        reset=false)
    
    reset && (xc = copy(booster.pos))

    for i in indices
        move(booster,x[:,i]; additive=false)
        updateHist!(booster,hist,freqs,objFunction; force=true)
    end

    reset && move(booster,xc; additive=false)

    return (a->a.objvalue).(hist[1:length(indices)])
end

const DefaultSimplexSampler = Callback(getSimplexObj)
