
using LinearAlgebra
using FFTW: fft, fft!, ifft, ifft!, fftshift, fftshift!, ifftshift, ifftshift!
using SpecialFunctions, FunctionZeros
using OffsetArrays: OffsetArray, OffsetMatrix, Origin, no_offset_view
using Plots

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

    @assert all(dz .> 0) "All propagated distances dz must be positive."

    eps = complex(eps)

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
    
    @assert all(dz .> 0) "All propagated distances dz must be positive."

    eps = complex(eps)

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
m = Modes(3,0,coords);


E0 = ones(C64,axes(coords.R));
E0 .*= coords.diskmaskin;


ax = axionModes(coords,m)
E = modes2field(ax,m)


coeffs = field2modes(E0,m)



freqs = collect(range(22e9,22.05e9,11))
d = collect(range(1e-3,10e-3,10))
@time p = propagationMatrix(d,freqs,1.0,m,coords);


include("spline.jl")









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


