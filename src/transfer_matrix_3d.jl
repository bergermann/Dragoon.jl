
using LinearAlgebra
using FFTW: fft, fft!, ifft, ifft!, fftshift, fftshift!, ifftshift, ifftshift!
using SpecialFunctions, FunctionZeros
using OffsetArrays: OffsetArray, OffsetMatrix, Origin
using Plots

const c0 = 299792458.


function kSpace(X)
    @assert X == -reverse(X) "Coordinates need to be symmetric around 0."

    k_max = π*length(X)/maximum(X)

    return ifftshift(range(-k_max,k_max,length(X)))
end

struct Coordinates
    X::Vector{Float64}
    kX::Vector{Float64}
    R::Matrix{Float64}
    kR::Matrix{Float64}
    Φ::Matrix{Float64}

    diskR::Float64
    diskmaskin::BitMatrix
    diskmaskout::BitMatrix

    function Coordinates(X::AbstractVector=-0.5:0.01:0.5; diskR::Real=0.15)
        @assert X == -reverse(X) "Coordinates must be symmetrical around 0."

        kX = kSpace(X)
        R  = [sqrt(x^2+y^2) for x in X, y in X]
        m  = R .<= diskR
        
        new(X,kX,R,
            [kx^2+ky^2 for kx in kX, ky in kX],
            [atan(y,x) for  x in  X,  y in  X],
            diskR,
            m,.!m
        )
    end

    function Coordinates(xsize::Real,dx::Real; diskR::Real=0.15)
        @assert xsize*dx > 0 "Inputs must be larger than 0."
    
        nx = ceil(xsize/2dx);
        X = -nx*dx:dx:nx*dx
        kX = kSpace(X)

        R  = [sqrt(x^2+y^2) for x in X, y in X]
        m  = R .<= diskR
    
        new(X,kX,R,
            [kx^2+ky^2 for kx in kX, ky in kX],
            [atan(y,x) for  x in  X,  y in  X],
            diskR,
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

    function Modes(M,L,coords)
        @assert M > 0 "m needs to be larger than 0."

        L_ = 2L+1
        patterns = Origin(1,-L,1,1,1)(
                zeros(ComplexF64,M,L_,length(coords.X),length(coords.Y),1))
        kt = Origin(1,-L)(zeros(ComplexF64,M,L_))

        display(typeof(patterns))
        display(typeof(kt))
        
        for m in 1:M, l in -L:L
            kt[m,l], patterns[m,l,:,:,:] = mode(m,l,coords)
        end

        return Modes(M,L,patterns,kt)
    end
end

function setMasks(coords::Coordinates,diskR::Real)
    coords.diskR = diskR
    coords.diskmaskin = coords.R .<= diskR
    coords.diskmaskout = .!coords.maskin

    return
end

function setMasks(coords::Coordinates)
    coords.diskmaskin = coords.R .<= coords.diskR
    coords.diskmaskout = .!coords.maskin

    return
end

function mode(m::Integer,l::Integer,coords::Coordinates)
    kr = besselj_zero(l,m)/diskR

    pattern = @. besselj(l,kr*coords.R)*cis(-l*coords.Φ)
    pattern[coords.R .> coords.diskR] .= 0.
    pattern ./= sqrt(sum(abs2.(pattern)))
    pattern = reshape(pattern,(size(pattern)...,1))

    return kr, pattern
end



# function propagate!(E0::Matrix{ComplexF64},coords::Coordinates,dz::Real,k0::Number)
#     fft!(E0)
#     @. E0 *= cis(-coords.kR*dz/2k0)
#     ifft!(E0)

#     return
# end

function propagate!(E0::Matrix{ComplexF64},coords::Coordinates,dz::Real,k0::Number)
    fft!(E0)
    @. @views E0 *= cis(sqrt(k0^2-coords.kR)*dz)
    ifft!(E0)

    return
end

# function propagateL!(E0::Matrix{ComplexF64},coords::Coordinates,dz::Real,k0::Number)
#     fft!(E0)
#     @. E0 *= cis(-conj(sqrt(k0^2-coords.kR))*dz)
#     ifft!(E0)

#     return
# end


f = 20e9; ω = 2π*f; λ = c0/f
eps = complex(1)
k0 = 2π*f/c0*sqrt(eps)

coords = Coordinates(1,λ/2; diskR=0.15)


E0 = ones(ComplexF64,axes(coords.R));
E0 .*= coords.diskmaskin;
# EL = copy(E0);

heatmap(coords.X,coords.X,abs2.(E0))


propagate!(E0,coords,1e-2,k0)
# propagateL!(EL,coords,1e-2,k0)



heatmap(coords.X,coords.X,abs2.(E0))
# heatmap(coords.X,coords.X,abs2.(EL))

heatmap(coords.X,coords.X,angle.(E0))
# heatmap(coords.X,coords.X,angle.(EL))
# heatmap(coords.X,coords.X,-angle.(EL))





