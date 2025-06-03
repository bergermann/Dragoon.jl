
using LinearAlgebra
using FFTW: fft, fft!, ifft, ifft!, fftshift, fftshift!, ifftshift, ifftshift!
using SpecialFunctions, FunctionZeros
using OffsetArrays: OffsetArray, OffsetMatrix, Origin, no_offset_view
using StaticArrays
using Plots

include("spline.jl")
include("transfer_matrix.jl")

const OM = OffsetMatrix; const OA = OffsetArray; const O = Origin; const raw = no_offset_view
const C64 = ComplexF64

# const c0 = 299792458.


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
    modes::OA{C64,5,Array{C64,5}}
    kt::OM{C64,Matrix{C64}}
    id::Matrix{C64}
    zero::Matrix{C64}

    function Modes(M,L,modes,kt)                
        @assert M > 0 "m needs to be larger than 0."

        s = M*(2L+1)
        id = Matrix{C64}(I,s,s)
        z = zeros(C64,s,s)

        new(M,L,modes,kt,id,z)
    end

    function Modes(M,L,coords)
        @assert M > 0 "m needs to be larger than 0."

        L_ = 2L+1
        modes = O(1,-L,1,1,1)(
                zeros(C64,M,L_,length(coords.X),length(coords.X),1))
        kt = O(1,-L)(zeros(C64,M,L_))
        
        for m in 1:M, l in -L:L
            kt[m,l], modes[m,l,:,:,:] = mode(m,l,coords)
        end

        return Modes(M,L,modes,kt)
    end
end

fieldDims(modes::Modes) = size(modes.modes,5)

import Base.getindex, Base.setindex!, Base.size, Base.axes
Base.getindex(m::Modes,inds...) = getindex(m.modes,inds...)
Base.getindex(m::Modes,ind1,ind2) = getindex(m.modes,ind1,ind2,:,:)
Base.getindex(m::Modes,ind1,ind2,ind3) = getindex(m.modes,ind1,ind2,:,:,ind3)
Base.setindex!(m::Modes,x,inds...) = setindex!(m.modes,x,inds...)
Base.setindex(m::Modes,ind1,ind2) = setindex(m.modes,ind1,ind2,:,:)
Base.setindex(m::Modes,ind1,ind2,ind3) = setindex(m.modes,ind1,ind2,:,:,ind3)
Base.size(m::Modes) = size(m.modes)
Base.size(m::Modes,d::Integer) = size(m.modes,d)
Base.axes(m::Modes,d::Integer) = axes(m.modes,d) 

function mode(m::Integer,l::Integer,coords::Coordinates)
    kr = besselj_zero(l,m)/coords.diskR

    mode = @. besselj(l,kr*coords.R)*cis(-l*coords.Φ)
    @. mode *= coords.diskmaskin
    mode ./= sqrt(sum(abs2.(mode)))
    mode = reshape(mode,(size(mode)...,1))

    return kr, mode
end



# function propagate!(E0::Matrix{C64},coords::Coordinates,dz::Real,k0::Number)
#     fft!(E0)
#     @. @views E0 *= cis(sqrt(k0^2-coords.kR)*dz)
#     ifft!(E0)

#     return
# end

function propagate!(E0::Matrix{C64},coords::Coordinates,dz::Real,k0::Number)
    fft!(E0)
    @. E0 *= cis(-conj(sqrt(k0^2-coords.kR))*dz)
    ifft!(E0)

    return
end




function modeDecomp(E::Union{Matrix{C64},Array{C64,3}},modes::Modes;)
    @assert size(E,1) == size(E,2) == size(modes,3) "Grids of field and modes don't match."
    @assert size(E,3) == size(modes,5) "Dimensionality of field and modes doesn't match."

    N = sqrt(sum(abs2.(E)))

    coeffs = O(1,-modes.L)(zeros(C64,modes.M,2modes.L+1))
    for m in 1:modes.M, l in -modes.L:modes.L
        coeffs[m,l] = sum(@. conj(modes.modes[m,l,:,:,:])*E)/N
    end

    return coeffs
end

const modeDecomposition = const decomp = const field2modes = modeDecomp



function axionModes(coords::Coordinates,modes::Modes; velocity_x::Real=0,f::Real=20e9)
    d = fieldDims(modes)
    Ea = zeros(C64,length(coords.X),length(coords.X),d)
    Ea[:,:,1+Int(d==3)] .= 1; Ea .*= coords.diskmaskin

    # Inaccuracies of the emitted fields: BField and Velocity Effects
    if velocity_x != 0
        k_a = 2π*f/c0 # k = 2pi/lambda (c/f = lambda)

        for (i,x) in enumerate(coords.X)
            Ea[i,:,:] .*= cis(k_a*x*velocity_x)
        end
    end

    return modeDecomp(Ea,modes)
end

function modes2field(coeffs::OM{C64,Matrix{C64}},modes::Modes)
    @assert all(size(coeffs) .<= size(modes)[1:2]) "More coefficients than modes available."

    field = zeros(Complex{Float64},size(modes,3),size(modes,4),size(modes,5))

    for m in 1:modes.M, l in -modes.L:modes.L
         @. field += coeffs[m,l]*modes.modes[m,l,:,:,:]
    end

    return field
end



function showField(E::Array{C64}; kwargs...)
    for i in axes(E,3); display(heatmap(abs.(E[:,:,i]); right_margin=4Plots.mm,kwargs...)); end
end

function showField(E::Array{C64},coords::Coordinates; kwargs...)
    for i in axes(E,3); display(heatmap(coords.X,coords.X,abs.(E[:,:,i]); right_margin=4Plots.mm,kwargs...)); end
end

function showCross(E::Array{C64}; kwargs...)
    n = size(E,1); @assert isodd(n) "Grid for E needs odd edge lengths."; n2 = div(n+1,2)
    for i in axes(E,3); display(plot(abs.(E[n2,:,i]); kwargs...)); end
end

function showCross(E::Array{C64},coords::Coordinates; kwargs...)
    n = size(E,1); @assert isodd(n) "Grid for E needs odd edge lengths."; n2 = div(n+1,2)
    for i in axes(E,3); display(plot(coords.X,abs.(E[n2,:,i]); kwargs...)); end
end



function propMatFreeSpace(freqs::AbstractVector{<:Real},distances::AbstractVector{<:Real},
        eps::Number,modes::Modes,coords::Coordinates)

    P = O(1,1, 1,-modes.L, 1,-modes.L)(
        zeros(C64,length(freqs),length(distances), modes.M,2modes.L+1, modes.M,2modes.L+1))

    for i in eachindex(freqs)
        k0 = 2π*freqs[i]/c0*sqrt(eps)

        for j in eachindex(distances)
            for m in 1:modes.M, l in -modes.L:modes.L, k in axes(modes,5)
                mode_ = copy(raw(modes[m,l,k]))
                propagate!(mode_,coords,distances[j],k0)
                coeffs_ = modeDecomp(mode_,modes)
                P[i,j,m,l,:,:] .+= coeffs_
            end
        end
    end

    return P
end

function propMatWaveGuide(freqs::AbstractVector{<:Real},distances::AbstractVector{<:Real},
        eps::Number,modes::Modes,coords::Coordinates)

    P = O(1,1, 1,-modes.L, 1,-modes.L)(
        zeros(C64,length(freqs),length(distances), modes.M,2modes.L+1, modes.M,2modes.L+1))

    for i in eachindex(freqs)
        k0 = 2π*freqs[i]/c0*sqrt(eps)

        for j in eachindex(distances)
            for m in 1:modes.M, l in -modes.L:modes.L#, k in axes(modes,5)
                # add disk tilts and surface here
                P[i,j,m,l,m,l] .= cis(-k0*distances[i])
            end
        end
    end

    return P
end


function propagationMatrix(freqs::AbstractVector{<:Real},distances::AbstractVector{<:Real},
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
    freqs::Vector{Float64}
    thickness::Float64
    nd::ComplexF64

    M::Int
    L::Int

    PS::OffsetArray{Spline{ComplexF64},5,Array{Spline{ComplexF64},5}}

    # work matrices for transfer_matrix algorithm
    Gd::SMatrix{2,2,ComplexF64}
    Gv::SMatrix{2,2,ComplexF64}
    G0::SMatrix{2,2,ComplexF64}
    #  T::MMatrix{2,2,ComplexF64}
     T::OffsetArray{ComplexF64,4,Array{ComplexF64,4}}

     S::SMatrix{2,2,ComplexF64}
    S0::SMatrix{2,2,ComplexF64}
    #  M::MMatrix{2,2,ComplexF64}
    MM::OffsetArray{ComplexF64,4,Array{ComplexF64,4}}

     W::MMatrix{2,2,ComplexF64}

    function GrandPropagationMatrix(freqs,distances,modes,coords; 
            eps::Real=24.0,tand::Real=0.0,thickness::Real=1e-3,nm::Real=1e15)

        p = propagationMatrix(freqs,distances,1.0,modes,coords);

        ps = O(1,1,-modes.L,1,-modes.L)(
        Array{Spline{C64}}(undef,length(freqs),modes.M,2modes.L+1,modes.M,2modes.L+1))

        for f in eachindex(freqs)
            for m in 1:modes.M, l in -modes.L:modes.L
                for m_ in 1:modes.M, l_ in -modes.L:modes.L
                    ps[f,m,l,m_,l_] = cSpline(d,p[f,:,m,l,m_,l_])
                end
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
        # T  = MMatrix{2,2,ComplexF64}(undef)
        T = O(1,1,1,-modes.L)(Array{ComplexF64}(undef,2,2,modes.M,2*modes.L+1))

        S  = SMatrix{2,2,ComplexF64}( A/2, 0.0im, 0.0im,  A/2)
        S0 = SMatrix{2,2,ComplexF64}(A0/2, 0.0im, 0.0im, A0/2)
        # M  = MMatrix{2,2,ComplexF64}(undef)
        MM = O(1,1,1,-modes.L)(Array{ComplexF64}(undef,2,2,modes.M,2*modes.L+1))

        W  = MMatrix{2,2,ComplexF64}(undef)

        new(freqs,thickness,nd,modes.M,modes.L,ps,Gd,Gv,G0,T,S,S0,MM,W)
    end
end

const GPM = GrandPropagationMatrix






# abstract type Space end
# abstract type Dist <: Space end
# abstract type Pos  <: Space end

function transfer_matrix_3d(::Type{Dist},distances::AbstractVector{<:Real},gpm::GPM,ax;)::OffsetArray{C64,4,Array{C64,4}}
    l = length(gpm.freqs)
    RB = O(1,1,1,-gpm.L)(Array{ComplexF64}(undef,l,2,gpm.M,2gpm.L+1))

    T_ = similar(gpm.T)
    W = gpm.W

    @inbounds @views for j in eachindex(gpm.freqs)
        f = gpm.freqs[j]

        pd1 = cispi(-2*f*gpm.nd*gpm.thickness/c0)
        pd2 = cispi(+2*f*gpm.nd*gpm.thickness/c0)

        for m in 1:gpm.M, l in -gpm.L:gpm.L
            copyto!(gpm.T[:,:,m,l],gpm.Gd)
            copyto!(gpm.MM[:,:,m,l],gpm.S); gpm.MM .*= ax[m,l]
        end

        # iterate in reverse order to sum up M in single sweep (thx david)
        for i in Iterators.reverse(eachindex(distances))
            for m in 1:gpm.M, l in -gpm.L:gpm.L
                gpm.T[:,1,m,l] .*= pd1
                gpm.T[:,2,m,l] .*= pd2 # T = Gd*Pd

                mul!(W,gpm.T[:,:,m,l],gpm.S); W .*= ax[m,l]; gpm.MM[:,:,m,l] .-= W         # M = Gd*Pd*S_-1
                mul!(W,gpm.T[:,:,m,l],gpm.Gv); copyto!(gpm.T[:,:,m,l],W)    # T *= Gd*Pd*Gv
            end

            T_ .= 0.0im
            for m in 1:gpm.M, l in -gpm.L:gpm.L                
                for m_ in 1:gpm.M, l_ in -gpm.L:gpm.L
                    s = spline(gpm.PS[j,m,l,m_,l_],distances[i])
                    # @assert abs(s) < 1.0 "|s| not smaller than 1!"
                    
                    T_[:,1,m,l] .+= gpm.T[:,1,m,l]*s
                    T_[:,2,m,l] .+= gpm.T[:,2,m,l]*conj(s)
                end
            end

            copyto!(gpm.T,T_)
               
            for m in 1:gpm.M, l in -gpm.L:gpm.L 
                if i > 1
                    mul!(W,gpm.T[:,:,m,l],gpm.S); W .*= ax[m,l]; gpm.MM[:,:,m,l] .+= W
                    mul!(W,gpm.T[:,:,m,l],gpm.Gd); copyto!(gpm.T[:,:,m,l],W)
                else
                    mul!(W,gpm.T[:,:,m,l],gpm.S0); W .*= ax[m,l]; gpm.MM[:,:,m,l] .+= W
                    mul!(W,gpm.T[:,:,m,l],gpm.G0); copyto!(gpm.T[:,:,m,l],W)
                end
            
                RB[j,1,m,l] = gpm.T[1,2,m,l]/gpm.T[2,2,m,l]
                RB[j,2,m,l] = gpm.MM[1,1,m,l]+gpm.MM[1,2,m,l]-
                    (gpm.MM[2,1,m,l]+gpm.MM[2,2,m,l])*gpm.T[1,2,m,l]/gpm.T[2,2,m,l]
            end
        end

        # M .= S
        # T .= 1.0+0.0im; T[1,1] += nd; T[2,2] += nd; T[2,1] -= nd; T[1,2] -= nd; T .*= 0.5
    end

    return RB
end


cispi(-2*f[1]*dists[1]/c0)
cispi( 2*f[1]*dists[1]/c0)
spline(gpm.PS[1,1,0,1,0],dists[1])
conj(spline(gpm.PS[1,1,0,1,0],dists[1]))


f = 22.025e9; ω = 2π*f; λ = c0/f
eps = complex(1)
k0 = 2π*f/c0*sqrt(eps)

coords = Coordinates(1,λ/2; diskR=0.15);
modes = Modes(3,0,coords);


# E0 = ones(C64,axes(coords.R));
# E0 .*= coords.diskmaskin;


# ax = axionModes(coords,modes)
# E = modes2field(ax,modes)


coeffs = field2modes(E0,modes)


# freqs = collect(range(21.5e9,22.5e9,101))
freqs = collect(range(21.98e9,22.26e9,50))
# @time p = propagationMatrix(freqs,collect(range(1e-3,10e-3,10)),1.0,modes,coords);
@time gpm = GPM(freqs,collect(range(1e-3,10e-3,10)),modes,coords; eps=24.0);



# dists = ones(20)*7.21e-3
dists = [1.00334,
        6.94754,
        7.1766,
        7.22788,
        7.19717,
        7.23776,
        7.07746,
        7.57173,
        7.08019,
        7.24657,
        7.21708,
        7.18317,
        7.13025,
        7.2198,
        7.45585,
        7.39873,
        7.15403,
        7.14252,
        6.83105,
        7.42282,]*1e-3

m_ = 3; l_ = 0;
@time RB = transfer_matrix_3d(Dist,dists,gpm,ax;);
@time RB_ = transfer_matrix(Dist,freqs,dists;); R_ = RB_[:,1]; B_ = RB_[:,2];

R = RB[:,1,m_,l_]; B = RB[:,2,m_,l_];
# display(plot(freqs/1e9,abs2.(B_),title="Boost 1d"))
display(plot(freqs/1e9,abs2.(B[:]*ax[1]),title="Boost 3d, m=$m_, l=$l_"))

# display(plot(freqs/1e9,abs.(R_),title="Ref 1d"))
# display(plot(freqs/1e9,abs.(R[:]),title="Ref 3d"))


gpm.PS[1,1,0,1,0]


# p = propagationMatrix(d,freqs,24,modes,coords)

# s = cSpline(d,p[100,:,1,0,1,0])
# spline(s,1.21e-3)


