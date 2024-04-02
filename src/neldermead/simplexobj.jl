###     functions to get objective values for every vertex

export getSimplexObj
export DefaultSimplexSampler

"""
    getSimplexObj(x::Matrix{Float64},
        booster::Booster,
        hist::Array{State},
        freqs::Array{Float64},
        objFunction::Callback,
        ()::Tuple{};
        reset=false)

Iteratively sample all vertices in ascending order, write result to x. Return to starting
position if `reset`.
"""
function getSimplexObj(x::Matrix{Float64},
                        booster::Booster,
                        hist::Array{State},
                        freqs::Array{Float64},
                        objFunction::Callback,
                        ()::Tuple;
                        reset=false)
    
    reset && (xc = copy(booster.pos))

    for i in axes(x,2)
        move(booster,x[:,i]; additive=false)
        updateHist!(booster,hist,freqs,objFunction; force=true)
    end

    reset && move(booster,xc; additive=false)

    return (a->a.objvalue).(hist[axes(x,2)])
end

"""
    getSimplexObj(x::Matrix{Float64},
        booster::Booster,
        hist::Array{State},
        freqs::Array{Float64},
        objFunction::Callback,
        ()::Tuple{};
        reset=false)

Iteratively sample all vertices by given indices in ascending order, write result to x.
Return to starting position if `reset`.
"""
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

    return reverse((a->a.objvalue).(hist[1:length(indices)]))
end

"""
    DefaultSimplexSampler

Callback for Nelder-Mead simplex sampler option by iteration for [`nelderMead`](@ref).
See [`getSimplexObj`](@ref).
"""
DefaultSimplexSampler = Callback(getSimplexObj)
