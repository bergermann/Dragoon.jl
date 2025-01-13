
using LinearAlgebra
using FFTW: fft, fft!, ifft, ifft!, fftshift, fftshift!, ifftshift, ifftshift!
using SpecialFunctions, FunctionZeros
using OffsetArrays: OffsetArray, OffsetMatrix, Origin
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

    diskmaskin::BitMatrix
    diskmaskout::BitMatrix

    function Coordinates(X::AbstractVector=-0.5:0.01:0.5,Y::AbstractVector=-0.5:0.01:0.5;
            diskR::Real=0.)
            
        kX = get_kspace_coords(X)
        kY = get_kspace_coords(Y)
        R  = [sqrt(x^2+y^2) for x in X, y in Y]
        m  = R .<= diskR
        
        new(X,Y,kX,kY,R,
            [kx^2+ky^2 for kx in kX, ky in kY],
            [atan(y,x) for  x in  X,  y in  Y],
            m,.!m
        )
    end

    function Coordinates(xsize::Real,dx::Real,ysize=nothing,dy=nothing; diskR::Real=0.)
        if isnothing(ysize); ysize=xsize; end
        if isnothing(dy);       dy=dx;    end

        @assert xsize*dx*ysize*dy > 0 "All inputs must be larger than 0."
    
        nx = ceil(xsize/2dx); X = -nx*dx:dx:nx*dx
        ny = ceil(ysize/2dy); Y = -ny*dy:dy:ny*dy

        kX = get_kspace_coords(X)
        kY = get_kspace_coords(Y)

        R  = [sqrt(x^2+y^2) for x in X, y in Y]
        m  = R .<= diskR
    
        new(X,Y,kX,kY,R,
        [kx^2+ky^2 for kx in kX, ky in kY],
        [atan(y,x) for  x in  X,  y in  Y],
            m,.!m
        )
    end
end



mutable struct Modes
    M::Int64
    L::Int64
    patterns::OffsetArray{ComplexF64,5,Array{ComplexF64,5}}
    kt::OffsetMatrix{ComplexF64,Matrix{ComplexF64}}
    id::Matrix{ComplexF64}
    zero::Matrix{ComplexF64}

    function Modes(M,L,patterns,kt)                
        @assert M > 0 "m needs to be larger than 0."

        s = M*(2L+1)
        id = Matrix{ComplexF64}(I,s,s)
        z = zeros(ComplexF64,s,s)

        new(M,L,patterns,kt,id,z)
    end

    function Modes(M,L,coords; diskR=0.15)
        @assert M > 0 "m needs to be larger than 0."

        L_ = 2L+1
        patterns = Origin(1,-L,1,1,1)(
                zeros(ComplexF64,M,L_,length(coords.X),length(coords.Y),1))
        kt = Origin(1,-L)(zeros(ComplexF64,M,L_))

        display(typeof(patterns))
        display(typeof(kt))
        
        for m in 1:M, l in -L:L
            kt[m,l], patterns[m,l,:,:,:] = mode(m,l,coords; diskR=diskR)
        end

        return Modes(M,L,patterns,kt)
    end
end

function setMasks(coords::Coordinates,diskR::Real)
    coords.diskmaskin = coords.R .<= diskR
    coords.diskmaskout = .!coords.maskin

    return
end

function mode(m::Integer,l::Integer,coords::Coordinates; diskR::Real=0.15)
    kr = besselj_zero(l,m)/diskR

    pattern = @. besselj(l,kr*coords.R)*cis(-l*coords.Φ)
    pattern[coords.R .> diskR] .= 0.
    pattern ./= sqrt(sum(abs2.(pattern)))
    pattern = reshape(pattern,(size(pattern)...,1))

    return kr, pattern
end


coords = Coordinates(1,0.01; diskR=0.15);
m = Modes(1,0,coords; diskR=0.15);




function propagate!(E0::Matrix{ComplexF64},coords::Coordinates,dz::Real,k0::Number)
    @time fft!(E0)
    # @time E0 .= fftshift(E0)

    @time @. E0 *= cis(-sqrt(k0^2-coords.kR)*dz)

    # @time E0 .= ifftshift(E0)
    @time ifft!(E0)

    return
end



E0 = ones(ComplexF64,axes(coords.R));
E0 .*= coords.diskmaskin; E0_ = copy(E0)

heatmap(abs2.(E0))

const c0 = 299792458.
f = 20e9
eps = complex(1)
k0 = 2π*f/c0*sqrt(eps)

propagate!(E0,coords,1e-1,k0)


heatmap(abs2.(E0))


