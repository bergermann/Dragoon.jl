export initSimplexCoord, initSimplexAffine

#constructors of the initial simplex

function initSimplexCoord(x0::Array{Float64},x::Matrix{Float64},d::Float64)
    x[:,:] = repeat(x0,1,length(x0)+1)

    for i in 1:length(x0)
        x[i,i+1] += d
    end
end

function initSimplexAffine(x0::Array{Float64},p::Matrix{Float64},a::Float64,
                                                                    b::Float64)
    x[:,:]
end