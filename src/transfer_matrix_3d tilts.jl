using LinearAlgebra
using FFTW: fft, fft!, ifft, ifft!, fftshift, fftshift!, ifftshift, ifftshift!
using SpecialFunctions, FunctionZeros
# using OffsetArrays: OffsetArray, OffsetMatrix, Origin, no_offset_view
using StaticArrays
using Plots
using BenchmarkTools

# include("spline.jl")
include("transfer_matrix.jl")

# const OM = OffsetMatrix; const OA = OffsetArray; const O = Origin; const raw = no_offset_view

# const c0 = 299792458.


function kSpace(X)
    @assert X == -reverse(X) "Coordinates need to be symmetric around 0."

    k_max = π*length(X)/2maximum(X)

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
        @assert X==-reverse(X) && X[1]*X[end]<0 "Coordinates must be symmetrical around 0."

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

mutable struct Modes
    M::Int64
    L::Int64
    modes::Array{ComplexF64,4}
    kt::Vector{ComplexF64}
    id::Matrix{ComplexF64}
    zero::Matrix{ComplexF64}

    function Modes(M,L,modes,kt)                
        @assert M > 0 "m needs to be larger than 0."

        ML = M*(2L+1)
        id = Matrix{ComplexF64}(I,ML,ML)
        z = zeros(ComplexF64,ML,ML)

        new(M,L,modes,kt,id,z)
    end

    function Modes(coords,M,L)
        @assert M > 0 "m needs to be larger than 0."

        ML = M*(2L+1)
        modes = zeros(ComplexF64,length(coords.X),length(coords.X),1,ML)
        kt = zeros(ComplexF64,ML)
        
        for m in 1:M, l in -L:L
            ml = modeidx(m,l,L)
            kt[ml], modes[:,:,:,ml] = mode(coords,m,l)
        end

        return Modes(M,L,modes,kt)
    end
end

fieldDims(modes::Modes) = size(modes.modes,5)

function modeidx(m::Int,l::Int,L::Int)
    return (m-1)*(2L+1)+l+L+1
end

function modeidx(ml::Int,L::Int)
    return div(ml-1,(2L+1))+1, (ml-1)%(2L+1)+-L
end

import Base.getindex, Base.setindex!, Base.size, Base.axes
Base.getindex(m::Modes,inds...) = getindex(m.modes,inds...)
Base.getindex(m::Modes,ind1,ind2) = getindex(m.modes,:,:,:,modeidx(ind1,ind2,m.L))
Base.getindex(m::Modes,ind1,ind2,ind3) = getindex(m.modes,:,:,ind1,modeidx(ind2,ind3,m.L))
Base.setindex!(m::Modes,x,inds...) = setindex!(m.modes,x,inds...)
Base.setindex(m::Modes,ind1,ind2) = setindex(m.modes,:,:,:,modeidx(ind1,ind2,m.L))
Base.setindex(m::Modes,ind1,ind2,ind3) = setindex(m.modes,:,:,ind1,modeidx(ind2,ind3,m.L))
Base.size(m::Modes) = size(m.modes)
Base.size(m::Modes,d::Integer) = size(m.modes,d)
Base.axes(m::Modes,d::Integer) = axes(m.modes,d) 

function mode(coords::Coordinates,m::Integer,l::Integer)
    kr = besselj_zero(l,m)/coords.diskR

    mode = @. besselj(l,kr*coords.R)*cis(-l*coords.Φ)
    @. mode *= coords.diskmaskin
    mode ./= sqrt(sum(abs2.(mode)))
    mode = reshape(mode,(size(mode)...,1))

    return kr, mode
end


function propagate!(E0::Matrix{C64},k0::Number,coords::Coordinates,dz::Real)
    fft!(E0)
    @. E0 *= cis(-conj(sqrt(k0^2-coords.kR))*dz)
    ifft!(E0)

    return
end

function propagate!(E0::Matrix{C64},k0::Number,coords::Coordinates,dz::Real,
        tiltx::Real,tilty::Real)

        
    fft!(E0)
    @. E0 *= cis(-conj(sqrt(k0^2-coords.kR))*dz)
    ifft!(E0)

    return
end



function modeDecomp(E::Union{Matrix{ComplexF64},Array{ComplexF64,3}},modes::Modes;)
    @assert size(E,1) == size(E,2) == size(modes,1) "Grids of field and modes don't match."
    @assert size(E,3) == size(modes,3) "Dimensionality of field and modes doesn't match."

    ML = modes.M*(2modes.L+1)

    N = sqrt(sum(abs2.(E)))
    coeffs = zeros(ComplexF64,ML)

    for ml in 1:ML; coeffs[ml] = sum(@. conj(modes.modes[:,:,:,ml])*E)/N; end

    return coeffs
end

const modeDecomposition = const decomp = const field2modes = modeDecomp



function axionModes(coords::Coordinates,modes::Modes; velocity_x::Real=0,f::Real=20e9)
    d = fieldDims(modes)
    Ea = zeros(ComplexF64,length(coords.X),length(coords.X),d)
    Ea[:,:,1+Int(d==3)] .= 1; Ea .*= coords.diskmaskin

    # inaccuracies of the emitted fields: B-field and velocity effects
    if velocity_x != 0
        k_a = 2π*f/c0 # k = 2pi/lambda (c/f = lambda)

        for (i,x) in enumerate(coords.X); Ea[i,:,:] .*= cis(k_a*x*velocity_x); end
    end

    return modeDecomp(Ea,modes)
end

function modes2field(coeffs::ComplexF64,modes::Modes)
    field = zeros(Complex{Float64},size(modes,1),size(modes,2),size(modes,3))

    for ml in 1:modes.M*(2modes.L+1); @. field += coeffs[ml]*modes.modes[:,:,:,ml]; end

    return field
end



function showField(E::Array{ComplexF64}; kwargs...)
    for i in axes(E,3)
        display(heatmap(abs.(E[:,:,i]); right_margin=4Plots.mm,kwargs...))
    end

    return
end

function showField(E::Array{ComplexF64},coords::Coordinates; kwargs...)
    for i in axes(E,3)
        display(heatmap(coords.X,coords.X,abs.(E[:,:,i]); right_margin=4Plots.mm,kwargs...))
    end

    return
end

function showCross(E::Array{ComplexF64}; kwargs...)
    n = size(E,1); @assert isodd(n) "Grid for E needs odd edge lengths."; n2 = div(n+1,2)
    for i in axes(E,3); display(plot(abs.(E[n2,:,i]); kwargs...)); end

    return
end

function showCross(E::Array{ComplexF64},coords::Coordinates; kwargs...)
    n = size(E,1); @assert isodd(n) "Grid for E needs odd edge lengths."; n2 = div(n+1,2)
    for i in axes(E,3); display(plot(coords.X,abs.(E[n2,:,i]); kwargs...)); end

    return
end



function propMatFreeSpace(freqs::Union{Real,AbstractVector{<:Real}},distances::AbstractVector{<:Real},
        eps::Number,modes::Modes,coords::Coordinates)

    ML = modes.M*(2modes.L+1)
    P = Array{ComplexF64}(undef,ML,ML,length(distances),length(freqs))

    for j in eachindex(freqs)
        k0 = 2π*freqs[j]/c0*sqrt(eps)

        for i in eachindex(distances)
            for ml in 1:ML#, k in axes(modes,3)
                mode = copy(modes[:,:,1,ml])            # mode = copy(modes[:,:,k,ml])
                propagate!(mode,k0,coords,distances[i]) # CHECK: kr here?
                coeffs = modeDecomp(mode,modes)
                @views copyto!(P[:,ml,i,j],coeffs)
            end
        end
    end

    return P
end

function propMatWaveGuide(freqs::AbstractVector{<:Real},distances::AbstractVector{<:Real},
        eps::Number,modes::Modes,coords::Coordinates)

    ML = modes.M*(2*modes.L+1)
    P = Array{ComplexF64}(undef,ML,ML,length(distances),length(freqs))

    for j in eachindex(freqs)
        k0 = 2π*freqs[i]/c0*sqrt(eps)

        for i in eachindex(distances)
            for ml in 1:ML
                P[ml,ml,i,j] .= cis(-k0*distances[i])
            end
        end
    end

    return P
end


function propagationMatrix(freqs::Union{Real,AbstractVector{<:Real}},distances::AbstractVector{<:Real},
        eps::Number,modes::Modes,coords::Coordinates; waveguide::Bool=false)

    @assert all(distances .> 0) "All propagated distances dz must be positive."
    @assert all(freqs .> 0) "All frequencies freqs must be positive."

    eps = complex(eps)

    if waveguide
        return propMatWaveGuide(freqs,distances,eps,modes,coords)
    else
        return propMatFreeSpace(freqs,distances,eps,modes,coords)
    end
end

const propMatrix = const propMat = const prop = propagationMatrix



mutable struct GrandPropagationMatrix
    freqs::Union{Real,Vector{Float64}}
    thickness::Float64
    nd::ComplexF64

    M::Int
    L::Int
    ML::Int

    PS::Array{Spline{ComplexF64},3}

    # work matrices for transfer_matrix algorithm
    Gd::SMatrix{2,2,ComplexF64}
    Gv::SMatrix{2,2,ComplexF64}
    G0::SMatrix{2,2,ComplexF64}
    
     S::SMatrix{2,2,ComplexF64}
    S0::SMatrix{2,2,ComplexF64}
    
     T::Vector{MMatrix{2,2,ComplexF64}}
    MM::Vector{MMatrix{2,2,ComplexF64}}

     W::MMatrix{2,2,ComplexF64}
    TW::Vector{MMatrix{2,2,ComplexF64}}

    function GrandPropagationMatrix(freqs,distances,modes,coords; 
            eps::Real=24.0,tand::Real=0.0,thickness::Real=1e-3,nm::Real=1e15)

        M = modes.M; L = 2modes.L+1; ML = M*L

        p = propagationMatrix(freqs,distances,1.0,modes,coords);

        ps = Array{Spline{ComplexF64}}(undef,ML,ML,length(freqs))

        for j in eachindex(freqs)
            for ml in 1:ML, ml_ in 1:ML
                ps[ml_,ml,j] = cSpline(distances,p[ml_,ml,:,j])
            end
        end

        ϵ  = eps*(1.0-1.0im*tand)
        nd = sqrt(ϵ); nm = complex(nm)
        ϵm = nm^2
        A  = 1-1/ϵ
        A0 = 1-1/ϵm

        Gd = SMatrix{2,2,ComplexF64}((1+nd)/2,   (1-nd)/2,   (1-nd)/2,   (1+nd)/2)
        Gv = SMatrix{2,2,ComplexF64}((nd+1)/2nd, (nd-1)/2nd, (nd-1)/2nd, (nd+1)/2nd)
        G0 = SMatrix{2,2,ComplexF64}((1+nm)/2,   (1-nm)/2,   (1-nm)/2,   (1+nm)/2)
        
        S  = SMatrix{2,2,ComplexF64}( A/2, 0.0im, 0.0im,  A/2)
        S0 = SMatrix{2,2,ComplexF64}(A0/2, 0.0im, 0.0im, A0/2)
        
        T  = [MMatrix{2,2,ComplexF64}(undef) for _ in 1:ML]
        MM = [MMatrix{2,2,ComplexF64}(undef) for _ in 1:ML]

        W  = MMatrix{2,2,ComplexF64}(undef)
        TW = [MMatrix{2,2,ComplexF64}(undef) for _ in 1:ML]

        new(freqs,thickness,nd,M,L,ML,ps,Gd,Gv,G0,S,S0,T,MM,W,TW)
    end
end

const GPM = GrandPropagationMatrix



# abstract type Space end
# abstract type Dist <: Space end
# abstract type Pos  <: Space end

function transfer_matrix_3d(::Type{Dist},distances::AbstractVector{<:Real},gpm::GPM,ax,f;)
    Gd = gpm.Gd; Gv = gpm.Gv; G0 = gpm.G0; S  = gpm.S; S0 = gpm.S0
    # T = gpm.T; TW = gpm.TW; MM  = gpm.MM; W = gpm.W

    PS = gpm.PS; M = gpm.M; L = 2modes.L+1; ML = M*L
    
    # RB = Array{ComplexF64}(undef,2,ML,length(gpm.freqs))
    RB = Array{ComplexF64}(undef,2,ML)
    
    # T  = [Matrix{ComplexF64}(undef,2,2) for _ in 1:ML]
    # MM = [Matrix{ComplexF64}(undef,2,2) for _ in 1:ML]
    T  = Array{ComplexF64}(undef,2,2,ML)
    MM = Array{ComplexF64}(undef,2,2,ML)

    W  = MMatrix{2,2,ComplexF64}(undef)
    # TW = [Matrix{ComplexF64}(undef,2,2) for _ in 1:ML]
    TW = Array{ComplexF64}(undef,2,2,ML)
    
    # @inbounds for j in eachindex(gpm.freqs)
    # @views @inbounds for j in eachindex(gpm.freqs)
        # f = gpm.freqs[j]
    # k = f*gpm.nd*gpm.thickness/c0
    pd1 = cispi(-2*f*gpm.nd*gpm.thickness/c0)
    pd2 = cispi(+2*f*gpm.nd*gpm.thickness/c0)
    # @time pd1 = cispi(-2*k)
    # @time pd2 = cispi(+2*k)

    # pd1 = exp(-2im*k*pi)
    # pd2 = exp(+2im*k*pi)

    for ml in eachindex(ax)
        copyto!(T[:,:,ml],Gd)
        copyto!(MM[:,:,ml],S)
    end
    
    MM .*= ax

    # iterate in reverse order to sum up MM in single sweep (thx david)
    @inbounds for i in Iterators.reverse(eachindex(distances))
        for ml in eachindex(ax)
            @. T[:,:,ml][:,1] *= pd1
            @. T[:,:,ml][:,2] *= pd2                     # T = Gd*Pd
            
            mul!(MM[:,:,ml],T[:,:,ml],S,-ax[ml],1.)              # MM = Gd*Pd*S_-1
            mul!(W,T[:,:,ml],Gv); copyto!(T[:,:,ml],W)            # T *= Gd*Pd*Gv

        end
        
        TW .= 0.0im

        for ml in eachindex(ax)
            for ml_ in eachindex(ax)
                s = spline(PS[ml,ml_,1],distances[i])

                @views @. TW[:,1,ml] += T[:,1,ml_]*s
                @views @. TW[:,2,ml] += T[:,2,ml_]*conj(s)
            end

            # copyto!(T[:,:,ml],TW[:,:,ml])
        # end
            
        # for ml in eachindex(ax)
            if i > 1
                mul!(MM[:,:,ml],T[:,:,ml],S,ax[ml],1.)
                mul!(W,T[:,:,ml],Gd); copyto!(T[:,:,ml],W)
            else
                mul!(MM[:,:,ml],T[:,:,ml],S0,ax[ml],1.)
                mul!(W,T[:,:,ml],G0); copyto!(T[:,:,ml],W)
            end
        
            # RB[1,ml,j] = T[:,:,ml][1,2]/T[:,:,ml][2,2]
            # RB[2,ml,j] = MM[:,:,ml][1,1]+MM[:,:,ml][1,2]-
            #     (MM[:,:,ml][2,1]+MM[:,:,ml][2,2])*T[:,:,ml][1,2]/T[:,:,ml][2,2]
        
            RB[1,ml] = T[:,:,ml][1,2]/T[:,:,ml][2,2]
            RB[2,ml] = MM[:,:,ml][1,1]+MM[:,:,ml][1,2]-
                (MM[:,:,ml][2,1]+MM[:,:,ml][2,2])*T[:,:,ml][1,2]/T[:,:,ml][2,2]
        end
    end
    # end

    return RB
end

@btime RB = transfer_matrix_3d(Dist,$dists,$gpm,$ax,$freqs;);

test();

function test()
    m = 1; l = 0
    
    freqs = 22.0e9; dists = ones(1)*7.21e-3
    coords = Coordinates(1,0.02; diskR=0.15);
    
    modes = Modes(coords,m,l); ax = axionModes(coords,modes)
    
    gpm = GPM(freqs,collect(range(1e-3,10e-3,10)),modes,coords; eps=24.0)
    @btime RB = transfer_matrix_3d(Dist,$dists,$gpm,$ax,$freqs;);
    
    return
end





f = 22.025e9; ω = 2π*f; λ = c0/f
eps = complex(1)
k0 = 2π*f/c0*sqrt(eps)

# coords = Coordinates(1,λ/2; diskR=0.15);
coords = Coordinates(1,0.02; diskR=0.15);

m = 1; l = 0
modes = Modes(coords,m,l);


ax = axionModes(coords,modes)

# freqs = collect(range(21.98e9,22.26e9,100));
freqs = 22.0e9

# @time p = propagationMatrix(freqs,collect(range(1e-3,10e-3,10)),1.0,modes,coords);
# @time gpm = GPM(freqs,collect(range(1e-3,10e-3,11)),modes,coords; eps=24.0);
@time gpm = GPM(freqs,collect(range(1e-3,10e-3,10)),modes,coords; eps=24.0);


dists = ones(1)*7.21e-3


# @time RB = transfer_matrix_3d(Dist,dists,gpm,ax;);
# RB = transfer_matrix_3d(Dist,dists,gpm,ax;);






B = [abs2.(RB[:,2,i,1]) for i in 1:m]
display(plot(freqs/1e9,B,title="Boost 3d, m_max = $m, l_max = $l",label=["m=1" "m=2" "m=3"]))


