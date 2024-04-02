###     constructors of the initial simplex

export initSimplexCoord, initSimplexAffine, initSimplexRegular
export InitSimplexCoord, InitSimplexAffine, InitSimplexRegular

"""
    initSimplexCoord(x0::Array{Float64},(d,)::Tuple{Float64})

Create initial simplex by moving `d` from `x0` along every coordinate axis.
"""
function initSimplexCoord(x0::Array{Float64},(d,)::Tuple{Float64})
    x = repeat(x0,1,length(x0)+1)

    for i in 1:length(x0)
        x[i,i+1] += d
    end

    return x
end


"""
    InitSimplexCoord(d)

Callback for Nelder-Mead simplex initialization along coordinate axes in [`nelderMead`](@ref).
See [`initSimplexCoord`](@ref).
"""
InitSimplexCoord(d) = Callback(initSimplexCoord,(d,))



"""
    initSimplexAffine(x0::Array{Float64},(a,b)::Tuple{Float64,Float64})

Create guaranteed affine initial simplex following [`this`](https://julianlsolvers.\
github.io/Optim.jl/v0.9.3/algo/nelder_mead/). 
"""
function initSimplexAffine(x0::Array{Float64},(a,b)::Tuple{Float64,Float64})
    x = repeat(x0,1,length(x0)+1)
    
    for i in eachindex(x0)
        x[i,i+1] *= b
        x[i,i+1] += a
    end

    return x
end

"""
InitSimplexAffine(a,b)

Callback for Nelder-Mead affine simplex initialization in [`nelderMead`](@ref).
See [`initSimplexAffine`](@ref).
"""
InitSimplexAffine(a,b) = Callback(initSimplexAffine,(a,b))


"""
    inCircleRadius(ndims::Int)

Return incircle radius of `n dimensional` regular polyhedron.
"""
function inCircleRadius(ndims::Int)
    return 1/sqrt(2*ndims*(ndims+1))
end

"""
    circumCircleRadius(ndims::Int)

Return circumcircle radius of `n dimensional` regular polyhedron.
"""
function circumCircleRadius(ndims::Int)
    return sqrt(ndims/(2*(ndims+1)))
end

"""
    initSimplexRegular(x0::Array{Float64},(d,)::Tuple{Float64,})

Create initial simplex as regular ndimensional polyhedron with edge length `d`.
"""
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

"""
    InitSimplexRegular(d)

Callback for regular Nelder-Mead simplex initialization in [`nelderMead`](@ref).
See [`initSimplexRegular`](@ref).
"""
InitSimplexRegular(d) = Callback(initSimplexRegular,(d,))



# simplex validation functions
"""
    simplexBaricenter(x::Matrix{Float64})

Return geometric baricenter/centroid of simplex.
"""
function simplexBaricenter(x::Matrix{Float64})
    return 1/size(x,2)*sum(x; dims=2)
end

"""
    simplexEdgeLengths(x::Matrix{Float64})

Return vector of all edge lengths of the simplex.
"""
function simplexEdgeLengths(x::Matrix{Float64})
    E = []

    for i in axes(x,2)
        for j in i+1:size(x,2)
            push!(E,pNorm(x[:,i]-x[:,j]))
        end
    end

    return E
end

"""
    simplexRadius(x::Matrix{Float64})

Return radii of all vertices from simplex baricenter.
"""
function simplexRadius(x::Matrix{Float64})
    x0 = simplexBaricenter(x)

    R = zeros(Float64,size(x,2))

    for i in axes(x,2)
        R[i] = pNorm(x0-x[:,i])
    end

    return R
end