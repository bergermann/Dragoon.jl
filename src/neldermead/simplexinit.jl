###     constructors of the initial simplex

export initSimplexCoord, initSimplexCoord!, initSimplexAffine
export InitSimplexCoord

# args = (d,)
function initSimplexCoord(x0::Array{Float64},args::Tuple{Float64})
    x = repeat(x0,1,length(x0)+1)

    for i in 1:length(x0)
        x[i,i+1] += args[1]
    end

    return x
end

const InitSimplexCoord(d) = Callback(initSimplexCoord,(d,))



# args = (p::Matrix{Float64},a::Float64,b::Float64,)
function initSimplexAffine(x0::Array{Float64},args::Tuple{})
    x = repeat(x0,1,length(x0)+1)
    x[:,:]

    return x
end
