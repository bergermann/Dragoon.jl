###     constructors of the initial simplex

export initSimplexCoord, initSimplexAffine
export InitSimplexCoord, InitSimplexAffine

# args = (d,)
function initSimplexCoord(x0::Array{Float64},args::Tuple{Float64})
    x = repeat(x0,1,length(x0)+1)

    for i in 1:length(x0)
        x[i,i+1] += args[1]
    end

    return x
end

const InitSimplexCoord(d) = Callback(initSimplexCoord,(d,))



# args = (a::Float64,b::Float64,)
function initSimplexAffine(x0::Array{Float64},args::Tuple{})
    x = repeat(x0,1,length(x0)+1)
    
    for i in eachindex(x0)
        x[i,i+1] *= b
        x[i,i+1] += a
    end

    return x
end

const InitSimplexAffine(a,b) = Callback(initSimplexCoord,(a,b))