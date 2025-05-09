
using LinearAlgebra
using FFTW: fft, fft!, ifft, ifft!, fftshift, fftshift!, ifftshift, ifftshift!
using SpecialFunctions, FunctionZeros
using OffsetArrays: OffsetArray, OffsetMatrix, Origin, no_offset_view
using StaticArrays
using Plots

include("spline.jl")

const OM = OffsetMatrix; const OA = OffsetArray; const O = Origin; const raw = no_offset_view
const C64 = ComplexF64

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



function propMatFreeSpace(dz::AbstractVector{<:Real},freqs::AbstractVector{<:Real},
        eps::Number,modes::Modes,coords::Coordinates)

    P = O(1,1, 1,-modes.L, 1,-modes.L)(
        zeros(C64,length(freqs),length(dz), modes.M,2modes.L+1, modes.M,2modes.L+1))

    for i in eachindex(freqs)
        k0 = 2π*freqs[i]/c0*sqrt(eps)

        for j in eachindex(dz)
            for m in 1:modes.M, l in -modes.L:modes.L, k in axes(modes,5)
                mode_ = copy(raw(modes[m,l,k]))
                propagate!(mode_,coords,dz[j],k0)
                coeffs_ = modeDecomp(mode_,modes)
                P[i,j,m,l,:,:] .+= coeffs_
            end
        end
    end

    return P
end

function propMatWaveGuide(dz::AbstractVector{<:Real},freqs::AbstractVector{<:Real},
        eps::Number,modes::Modes,coords::Coordinates)

    P = O(1,1, 1,-modes.L, 1,-modes.L)(
        zeros(C64,length(freqs),length(dz), modes.M,2modes.L+1, modes.M,2modes.L+1))

    for i in eachindex(freqs)
        k0 = 2π*freqs[i]/c0*sqrt(eps)

        for j in eachindex(dz)
            for m in 1:modes.M, l in -modes.L:modes.L#, k in axes(modes,5)
                # add disk tilts and surface here
                P[i,j,m,l,m,l] .= cis(-k0*dz[i])
            end
        end
    end

    return P
end


function propagationMatrix(dz::AbstractVector{<:Real},freqs::AbstractVector{<:Real},
        eps::Number,modes::Modes,coords::Coordinates; waveguide::Bool=false)

    @assert all(dz .> 0) "All propagated distances dz must be positive."
    @assert all(freqs .> 0) "All frequencies freqs must be positive."

    eps = complex(eps)

    if waveguide
        return propMatWaveGuide(dz,freqs,eps,modes,coords)
    else
        return propMatFreeSpace(dz,freqs,eps,modes,coords)
    end
end

const propMatrix = const propMat = const prop = propagationMatrix



f = 20e9; ω = 2π*f; λ = c0/f
eps = complex(1)
k0 = 2π*f/c0*sqrt(eps)

coords = Coordinates(1,λ/2; diskR=0.15);
modes = Modes(3,0,coords);


E0 = ones(C64,axes(coords.R));
E0 .*= coords.diskmaskin;


ax = axionModes(coords,modes)
E = modes2field(ax,modes)


coeffs = field2modes(E0,modes)



freqs = collect(range(22e9,22.05e9,11))
d = collect(range(1e-3,10e-3,10))
@time p = propagationMatrix(d,freqs,1.0,modes,coords);



# f_ = 1; m = 3; l = 0; m_ = 1; l_ = 0
# p_ = p[f_,:,m,l,m_,l_]
# s = cSpline(d,p_)

# xs = collect(1:0.01:10)*1e-3
# ys = spline(s,xs)
# plot(d/1e-3,[real.(p_),imag.(p_)]; seriestype=:scatter,label=permutedims(["re(data)","im(data)"]),
#     xlabel="d [mm]",title="F: $(freqs[f_]/1e9) GHz, m: $m, l: $l, m': $m_, l': $l_")
# plot!(xs/1e-3,[real.(ys),imag.(ys)]; label=permutedims(["re(spline)","im(spline)"]))



ps = O(1,1,-modes.L,1,-modes.L)(
        Array{Spline{C64}}(undef,length(freqs),modes.M,2modes.L+1,modes.M,2modes.L+1))

for f in eachindex(freqs)
    for m in 1:modes.M, l in -modes.L:modes.L
        for m_ in 1:modes.M, l_ in -modes.L:modes.L
            ps[f,m,l,m_,l_] = cSpline(d,p[f,:,m,l,m_,l_])
        end
    end
end

# for i in axes(a,1)
#     for m in 1:3, m_ in 1:3
#         if m == m_; continue; end
#         display(abs(a[i,m,m_]-a[i,m_,m]))
#         display(angle(a[i,m,m_])-angle(a[i,m_,m]))
#     end
# end


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
     T::OffsetArray{ComplexF64,4,Array{Complex64,4}}

     S::SMatrix{2,2,ComplexF64}
    S0::SMatrix{2,2,ComplexF64}
    #  M::MMatrix{2,2,ComplexF64}
     M::OffsetArray{ComplexF64,4,Array{Complex64,4}}

     W::MMatrix{2,2,ComplexF64}

    function GrandPropagationMatrix(freqs,distances,modes,coords; 
            eps::Real=24.0,tand::Real=0.0,thickness::Real=1e-3,nm::Real=1e15)

        p = propagationMatrix(distances,freqs,1.0,modes,coords);

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
        T = O(1,1,1,-modes.L)(Matrix{Complex64}(undef,2,2,modes.M,2*modes.L+1))

        S  = SMatrix{2,2,ComplexF64}( A/2, 0.0im, 0.0im,  A/2)
        S0 = SMatrix{2,2,ComplexF64}(A0/2, 0.0im, 0.0im, A0/2)
        # M  = MMatrix{2,2,ComplexF64}(undef)
        M = O(1,1,1,-modes.L)(Matrix{Complex64}(undef,2,2,modes.M,2*modes.L+1))

        W  = MMatrix{2,2,ComplexF64}(undef)

        new(freqs,thickness,nd,modes.M,modes.L,ps,Gd,Gv,G0,T,S,S0,M,W )
    end
end

const GPM = GrandPropagationMatrix




gpm = GPM(freqs,d,modes,coords)





abstract type Space end
abstract type Dist <: Space end
abstract type Pos  <: Space end

function transfer_matrix_3d(::Type{Dist},distances::AbstractVector{<:Real},gpm::GPM;)::Matrix{ComplexF64}
    l = length(gpm.freqs)
    RB = O(1,1,1,-gpm.L)(Array{ComplexF64}(undef,l,2,gpm.M,2gpm.L+1))



    T .= gpm.Gd
    M .= gpm.S

    @inbounds @views for j in eachindex(gpm.freqs)
        f = gpm.freqs[j]

        pd1 = cispi(-2*f*gpm.nd*gpm.thickness/c0)
        pd2 = cispi(+2*f*gpm.nd*gpm.thickness/c0)

        for m in 1:gpm.M, l in -gpm.L:gpm.L
            gpm.T[:,:,m,l] .= gpm.Gd
            gpm.M[:,:,m,l] .= gpm.S
        end

        # iterate in reverse order to sum up M in single sweep (thx david)
        for i in Iterators.reverse(eachindex(distances))
            for m in 1:gpm.M, l in -gpm.L:gpm.L
                T[:,1] .*= pd1
                T[:,2] .*= pd2 # T = Gd*Pd

                mul!(W,T,S); M .-= W    # M = Gd*Pd*S_-1

                mul!(W,T,Gv); T .= W    # T *= Gd*Pd*Gv

                T[:,1] .*= cispi(-2*f*distances[i]/c0)
                T[:,2] .*= cispi(+2*f*distances[i]/c0)   # T = Gd*Pd*Gv*Gd*S_-1

                if i > 1
                    mul!(W,T,S); M .+= W
                    mul!(W,T,Gd); T .= W
                else
                    mul!(W,T,S0); M .+= W
                    mul!(W,T,G0); T .= W
                end
            end
            
            RB[j] = T[1,2]/T[2,2]
            RB[l+j] = M[1,1]+M[1,2]-(M[2,1]+M[2,2])*T[1,2]/T[2,2]

            M .= S
            T .= 1.0+0.0im; T[1,1] += nd; T[2,2] += nd; T[2,1] -= nd; T[1,2] -= nd; T .*= 0.5
        end
    end

    return RB
end








function transformer(bdry::SetupBoundaries, coords::CoordinateSystem, modes::Modes; f=10.0e9, velocity_x=0, prop=propagator, propagation_matrices::Array{Array{Complex{T},2},1}=Array{Complex{Float64},2}[], diskR=0.15, emit=axion_induced_modes(coords,modes;B=nothing,velocity_x=velocity_x,diskR=diskR,f=f), reflect=nothing) where T<:Real
    bdry.eps[isnan.(bdry.eps)] .= 1e30

    transmissionfunction_complete = [modes.id modes.zeromatrix ; modes.zeromatrix modes.id ]
    lambda = wavelength(f)

    initial = emit

    axion_beam = Array{Complex{T}}(zeros((modes.M)*(2modes.L+1)))

    Nregions = length(bdry.eps)
    idx_reg(s) = Nregions-s+1

    for s in (Nregions-1):-1:1
        axion_beam .+= axion_contrib(transmissionfunction_complete, sqrt(bdry.eps[idx_reg(s+1)]), sqrt(bdry.eps[idx_reg(s)]), initial, modes)

        diffprop = (isempty(propagation_matrices) ?
                        propagation_matrix(bdry.distance[idx_reg(s)], diskR, bdry.eps[idx_reg(s)], bdry.relative_tilt_x[idx_reg(s)], bdry.relative_tilt_y[idx_reg(s)], bdry.relative_surfaces[idx_reg(s),:,:], lambda, coords, modes; prop=prop) :
                        propagation_matrices[idx_reg(s)])

        transmissionfunction_complete *= get_boundary_matrix(sqrt(bdry.eps[idx_reg(s)]), sqrt(bdry.eps[idx_reg(s+1)]), diffprop, modes)
    end

    boost =  - (transmissionfunction_complete[index(modes,2),index(modes,2)]) \ (axion_beam)

    if reflect === nothing
        return boost
    end

    refl = - transmissionfunction_complete[index(modes,2),index(modes,2)] \
           ((transmissionfunction_complete[index(modes,2),index(modes,1)]) * (reflect))
    return boost, refl
end



function propagation_matrix(dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem, modes::Modes; is_air=(real(eps)==1), onlydiagonal=false, prop=propagator)
    matching_matrix = Array{Complex{Float64}}(zeros(modes.M*(2modes.L+1),modes.M*(2modes.L+1)))


    k0 = 2pi/lambda*sqrt(eps)

    # Define the propagation function
    propfunc = nothing # initialize
    if is_air
        # In air use the propagator we get
        function propagate(x)
            return prop(copy(x), dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords)
        end
        propfunc = propagate
    else
        # In the disk the modes are eigenmodes, so we only have to apply the
        # inaccuracies and can apply the propagation later separately
        propfunc(efields) = efields.*[exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in coords.X, y in coords.Y].*exp.(-1im*k0*surface)
        # Applying exp element-wise to the surface is very important otherwise it is e^M with M matrix
    end

    # Calculate the mixing matrix
    for m_prime in 1:modes.M, l_prime in -modes.L:modes.L
        # Propagated fields of mode (m_prime, l_prime)
        for i in 1:e_field_dimensions(modes)
            # If 3D E-field, fields propagate separately. For interaction need to implement in propagator.
            propagated = propfunc(modes.mode_patterns[m_prime,l_prime+modes.L+1,:,:,i])

            for m in (onlydiagonal ? [m_prime] : 1:modes.M), l in (onlydiagonal ? [l_prime] : -modes.L:modes.L)
                # P_ml^m'l' = <E_ml | E^p_m'l' > = ∫ dA \bm{E}_ml^* ⋅ \bm{E}^p_m'l' = ∑_{j = x,y,z} ∫ dA ({E}_j)_ml^* ⋅ ({E}_j)^p_m'l'
                # 6.15 in Knirck
                matching_matrix[(m-1)*(2modes.L+1)+l+modes.L+1, (m_prime-1)*(2modes.L+1)+l_prime+modes.L+1] +=
                        sum( conj.(modes.mode_patterns[m,l+modes.L+1,:,:,i]) .* propagated ) #*dx*dy

                #v = 1-abs2.(matching_matrix[(m-1)*(2L+1)+l+L+1, (m_prime-1)*(2L+1)+l_prime+L+1])
            end
        end
    end

    if !is_air
        propagation_matrix = Array{Complex{Float64}}(zeros(modes.M*(2modes.L+1),modes.M*(2modes.L+1)))
        #The propagation within the disk is still missing
        for m in 1:modes.M, l in -modes.L:modes.L
            kz = sqrt(k0^2 - modes.mode_kt[m,l+modes.L+1]^2)
            propagation_matrix[(m-1)*(2modes.L+1)+l+modes.L+1, (m-1)*(2modes.L+1)+l+modes.L+1] = exp(-1im*kz*dz)
        end
        # It is important to note the multiplication from the left
        matching_matrix = propagation_matrix*matching_matrix
    end

    return matching_matrix

    # TODO: This only takes surface roughness at the end of the propagation into account, not at its
    # start. General problem in all the codes so far.
end


