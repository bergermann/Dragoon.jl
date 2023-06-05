###     constructors of the initial simplex

export initSimplexCoord, initSimplexCoord!, initSimplexAffine
export InitSimplexCoord


function initSimplexCoord(x0::Array{Float64},d::Real,args::Tuple{})
    x = repeat(x0,1,length(x0)+1)

    for i in 1:length(x0)
        x[i,i+1] += args[1]
    end

    return x
end

const InitSimplexCoord(d) = Callback(initSimplexCoord,(d,))

function initSimplexAffine(x0::Array{Float64},p::Matrix{Float64},a::Float64,
                                                    b::Float64,args::Tuple{})
    x[:,:]
end
