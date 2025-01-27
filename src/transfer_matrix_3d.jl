
using LinearAlgebra
using FFTW: fft, fft!, ifft, ifft!, fftshift, fftshift!, ifftshift, ifftshift!
using SpecialFunctions, FunctionZeros
using OffsetArrays: OffsetArray, OffsetMatrix, Origin
using Plots

const OM = OffsetMatrix; OA = OffsetArray; O = Origin

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
    modes::OA{ComplexF64,5,Array{ComplexF64,5}}
    kt::OM{ComplexF64,Matrix{ComplexF64}}
    id::Matrix{ComplexF64}
    zero::Matrix{ComplexF64}

    function Modes(M,L,modes,kt)                
        @assert M > 0 "m needs to be larger than 0."

        s = M*(2L+1)
        id = Matrix{ComplexF64}(I,s,s)
        z = zeros(ComplexF64,s,s)

        new(M,L,modes,kt,id,z)
    end

    function Modes(M,L,coords)
        @assert M > 0 "m needs to be larger than 0."

        L_ = 2L+1
        modes = O(1,-L,1,1,1)(
                zeros(ComplexF64,M,L_,length(coords.X),length(coords.X),1))
        kt = O(1,-L)(zeros(ComplexF64,M,L_))
        
        for m in 1:M, l in -L:L
            kt[m,l], modes[m,l,:,:,:] = mode(m,l,coords)
        end

        return Modes(M,L,modes,kt)
    end
end

fieldDims(modes::Modes) = size(modes.modes,5)

import Base.getindex, Base.setindex!
Base.getindex(m::Modes,inds...) = getindex(m.modes,inds...)
# Base.setindex(m::Modes,x,inds...) = setindex(m.modes,x,inds...)

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
    kr = besselj_zero(l,m)/coords.diskR

    mode = @. besselj(l,kr*coords.R)*cis(-l*coords.Φ)
    @. mode *= coords.diskmaskin
    mode ./= sqrt(sum(abs2.(mode)))
    mode = reshape(mode,(size(mode)...,1))

    return kr, mode
end


# function propagate!(E0::Matrix{ComplexF64},coords::Coordinates,dz::Real,k0::Number)
#     fft!(E0)
#     @. @views E0 *= cis(sqrt(k0^2-coords.kR)*dz)
#     ifft!(E0)

#     return
# end

function propagateL!(E0::Matrix{ComplexF64},coords::Coordinates,dz::Real,k0::Number)
    fft!(E0)
    @. E0 *= cis(-conj(sqrt(k0^2-coords.kR))*dz)
    ifft!(E0)

    return
end


f = 20e9; ω = 2π*f; λ = c0/f
eps = complex(1)
k0 = 2π*f/c0*sqrt(eps)

coords = Coordinates(1,λ/2; diskR=0.15);
m = Modes(3,3,coords);


E0 = ones(ComplexF64,axes(coords.R));
E0 .*= coords.diskmaskin;






function axionModes(coords::Coordinates,modes::Modes; B=nothing,velocity_x::Real=0,f::Real=20e9)
    if isnothing(B)
        d = fieldDims(modes)
        B = zeros(ComplexF64,length(coords.X),length(coords.X),d)
        B[:,:,1+Int(d==3)] .= 1
    end

    # Inaccuracies of the emitted fields: BField and Velocity Effects ###################
    if velocity_x != 0
        B = Array{Complex{Float64}}(B)
        k_a = 2π*f/c0 # k = 2pi/lambda (c/f = lambda)

        for (i,x) in enumerate(coords.X)
            B[i,:,:] .*= cis(k_a*x*velocity_x)
        end
    end

    # Only the axion-induced field on the disks matters:
    B .*= coords.diskmaskin

    # Note: in all the normalizations the dx factors are dropped, since they should drop out in the final
    #       integrals
    B ./= sqrt(sum(abs2.(B)))

    modes_ = O(1,-modes.L)(zeros(ComplexF64,modes.M,2modes.L+1))
    for m in 1:modes.M, l in -modes.L:modes.L
        modes_[m,l] = sum(@. B*conj(modes.modes[m,l,:,:,:]))
    end

    return modes_
end

function modes2field(mode_coeffs::OM{ComplexF64,Matrix{ComplexF64}},coords::Coordinates,modes::Modes)
    field = zeros(Complex{Float64},length(coords.X),length(coords.X),fieldDims(modes))

    for m in 1:modes.M, l in -modes.L:modes.L
         @. field += mode_coeffs[m,l]*modes.modes[m,l,:,:,:]
    end

    return field
end

field2modes(E::Array{<:ComplexF64},coords::Coordinates,modes::Modes) = axionModes(coords,modes; B=E)

ax = axionModes(coords,m)
E = modes2field(ax,coords,m)

heatmap(abs2.(E[:,:,1]))


coeffs = field2modes(E0,coords,m)
heatmap(abs2.(modes2field(coeffs,coords,m)[:,:,1]))

function showField(E::Array{ComplexF64})
    for i in axes(E,3)
        heatmap(abs2.(E[:,:,i]))
    end
end