
using LinearAlgebra

function eval_polynomial(x::Number,coeffs::AbstractVector)
    p = 0

    for i in eachindex(coeffs)
        p += coeffs[i]*x^(i-1)
    end

    return p
end

struct Spline
    order::Int
    knots::Vector{Float64}
    coeffs::Matrix{Float64}

    function Spline(order,knots,coeffs)
        @assert order+1 == size(coeffs,2)
        @assert length(knots) == size(coeffs,1)+1

        new(order,knots,coeffs)
    end
end

function cSpline(x,y; b=[1,0,0,1],c=[0,0])
    @assert length(x) == length(y) "x and y need same lengths."

    h = x[2]-x[1]
    n = length(x)-1

    A = zeros(n,4); A[:,1] .= y[1:end-1]

    du = ones(Float64,n-1); du[1] = b[2]
    d  = ones(Float64,n)*4; d[1]  = b[1]; d[end] = b[4]
    dl = ones(Float64,n-1); dl[end] = b[3]

    A2 = Tridiagonal(dl, d, du); A2_ = inv(A2)

    D = zeros(Float64,n); D[1] = c[1]-y[1]; D[end] = c[2]-y[end]

    for i in 2:n-1
        D[i-1] += y[i]
        D[i] -= 2y[i]
        D[i+1] += y[i]
    end

    D .*= 3/h

    A[:,3] = A2_*D
    @. A[1:end-1,4] = A[2:end,3]-A[1:end-1,3]
    @. A[:,2] = (y[2:end]-y[1:end-1])/h + A[:,3]*h + A[:,4]*h^2

    return Spline(3,x,A)
end

function cSpline(bounds,len,f; kwargs...)
    x = collect(range(bounds[1],bounds[2],len))
    y = f.(x)

    return cSpline(x,y; kwargs...)
end

function differentiate(spline::Spline)
    if spline.order == 0
        return Spline(0,spline.knots,zeros(size(spline.coeffs)))
    end

    coeffs = Matrix{undef,size(spline.coeffs,1),spline.order}

    for i in 1:spline.order
        coeffs[:,i] = i*spline.coeffs[:,i+1]
    end

    return Spline(spline.order-1,spline.knots,coeffs)
end

const diff = differentiate

function spline(spline::Spline,x::Real)
    @assert spline.knots[1] <= x <= spline.knots[end] "x needs to be within spline bounds."

    idx = min(findlast(y->y<=x,spline.knots),length(spline.knots)-1)

    return eval_polynomial(x-spline.knots[idx],spline.coeffs[idx,:])
end


x = collect(0:10)
y = rand(length(x))

S = cSpline(x,y)

spline(S,4)
y[5]

x_ = collect(0:0.01:10)
S_ = [spline(S,x_[i]) for i in eachindex(x_)]

scatter(x,y)
plot!(x_,S_)