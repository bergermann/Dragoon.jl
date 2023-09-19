###     constructors of the initial simplex

export initSimplexCoord, initSimplexAffine, initSimplexRegular
export InitSimplexCoord, InitSimplexAffine, InitSimplexRegular

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

const InitSimplexAffine(a,b) = Callback(initSimplexAffine,(a,b))



function inCircleRadius(ndims::Int)
    return 1/sqrt(2*ndims*(ndims+1))
end
    
function circumCircleRadius(ndims::Int)
    return sqrt(ndims/(2*(ndims+1)))
end

function initSimplexRegular(x0::Array{Float64},(d,)::Tuple{Float64,})
    x = zeros(length(x0),length(x0)+1)

    for i in eachindex(x0)
        x[i,1:i] .-= d*inCircleRadius(i)
        x[i,i+1] += d*circumCircleRadius(i)
    end

    for i in 1:length(x0)+1
        x[:,i] += x0
    end

    return x
end

const InitSimplexRegular(d) = Callback(initSimplexRegular,(d,))



# simplex validation functions
function simplexBaricenter(x::Matrix{Float64})
    return 1/size(x,2)*sum(x; dims=2)
end

function simplexEdgeLengths(x::Matrix{Float64})
    E = []

    for i in axes(x,2)
        for j in i+1:size(x,2)
            push!(E,pNorm(x[:,i]-x[:,j]))
        end
    end

    return E
end

function simplexRadius(x::Matrix{Float64})
    x0 = simplexBaricenter(x)

    R = zeros(Float64,size(x,2))

    for i in axes(x,2)
        R[i] = pNorm(x0-x[:,i])
    end

    return R
end