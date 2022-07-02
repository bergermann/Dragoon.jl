###     functions to get objective values for every vertex

export getSimplexObj

function getSimplexObj(x::Matrix{Float64},
                        indices::Vector{Int},
                        booster::Booster,
                        hist::Array{State},
                        freqs::Array{Float64},
                        objFunction::Tuple{Function,Vector};
                        reset=false)
    reset && (xc = copy(booster.pos))

    for i in indices
        moveCommand(booster,x[:,i]; additive=false)
        updateHist!(booster,hist,freqs,objFunction)
    end

    reset && moveCommand(booster,xc; additive=false)

    return (a->a.objvalue).(hist[1:length(indices)])
end
