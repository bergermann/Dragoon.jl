
using FFTW
using SpecialFunctions, FunctionZeros
using Plots



function get_kspace_coords(coords)
    min_k = π/maximum(coords)
    max_k = min_k*(length(coords)-1)/2
    coordsk = -max_k:min_k:max_k

    return coordsk
end

struct Coordinates
    X::Vector{Float64}
    Y::Vector{Float64}
    kX::Vector{Float64}
    kY::Vector{Float64}
    R::Matrix{Float64}
    kR::Matrix{Float64}
    Φ::Matrix{Float64}

    function Coordinates(X=-0.5:0.01:0.5,Y=-0.5:0.01:0.5)
        kX = get_kspace_coords(X)
        kY = get_kspace_coords(Y)
        
        new(X,Y,kX,kY,
            [  sqrt(x^2+y^2) for  x in  X,  y in  Y],
            [sqrt(kx^2+ky^2) for kx in kX, ky in kY],
            [      atan(x/y) for  x in  X,  y in  Y],
        )
    end

    function Coordinates(xsize,dx,ysize,dy)
        @assert xsize*dx*ysize*dy > 0 "All inputs must be larger than 0."
    
        nx = ceil(xsize/2dx); X = -nx*dx:dx:nx*dx
        ny = ceil(ysize/2dy); Y = -ny*dy:dy:ny*dy

        kX = get_kspace_coords(X)
        kY = get_kspace_coords(Y)
    
        new(X,Y,kX,kY,
            [  sqrt(x^2+y^2) for  x in  X,  y in  Y],
            [sqrt(kx^2+ky^2) for kx in kX, ky in kY],
            [      atan(x/y) for  x in  X,  y in  Y],
        )
    end
end



function mode(m,l,L,coords::Coordinates; diskR=0.15,pattern_input=nothing,kT_arr=nothing)    
    if isnothing(pattern_input)
        kr = besselj_zero(abs.(l),m)/diskR

        pattern = @. besselj(l,kr*coords.R)*cis(-l*coords.Φ)
        pattern[coords.R .> diskR] .= 0.; display(abs.(pattern))
        println(sum(abs2.(pattern)))
        pattern ./= sqrt(sum(abs2.(pattern))); display(abs.(pattern))
        pattern = reshape(pattern, (size(pattern)...,1)); display(abs.(pattern))
    else
        @assert !isnothing(kT_arr) "kTs need to be defined alongside pattern_input."
        
        pattern = pattern_input[m,l+L+1,:,:,:]
        kr = kT_arr[m,l+L+1]
    end

    return kr, pattern
end

coords = Coordinates(1,0.1,1,0.1)
mode(0,0,0,coords)



