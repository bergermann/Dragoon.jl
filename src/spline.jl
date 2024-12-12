
using StaticArrays

function eval_polynomial(x::Number,coeffs::AbstractVector)
    p = 0

    for i in eachindex(coeffs)
        p += coeffs[i]*x^i
    end

    return p
end

struct Spline
    order::Int
    knots::Vector{FLoat64}
    coeffs::Matrix{FLoat64}

    function Spline(order,knots,coeffs)
        @assert order+1 == size(coeffs,2)
        @assert length(order) == size(coeffs,1)+1

        new(order,knots,coeffs)
    end
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

function spline(spline::Spline,x::Float64)
    @assert spline.knots[1] <= x <= spline.knots[end] "x needs to be within spline bounds."

    idx = min(findlast(y->y<=x,spline.knots),length(spline.knots)-1)

    return eval_polynomial(x-spline.knots[idx],spline.coeffs[idx])
end


x = collect(0:10)
y = rand(length(x))
n = length(x)

coeffs = zeros(n-1,3)

S = Spline{Float64,3,n,n-1}(x,coeffs)