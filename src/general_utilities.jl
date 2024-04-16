###     convenient convenience functions for convenience

export spacings_stefan,
    boost1d, ref1d, getBoost1d, getRef1d,
    boost3d, getBoost3d,
    findpeak1d, findpeak3d,
    genFreqs, 
    pos2dist, dist2pos,
    pNorm, copy

"""
    const spacings_stefan

Good standard set of spacings by Stefan Knirck, see
[`BoostFractor repository`](https://github.com/mppmu/BoostFractor.jl).
"""
const spacings_stefan = [1.00334, 6.94754, 7.1766, 7.22788, 7.19717,
        7.23776, 7.07746, 7.57173, 7.08019, 7.24657,
        7.21708, 7.18317, 7.13025, 7.2198, 7.45585,
        7.39873, 7.15403, 7.14252, 6.83105, 7.42282]*1e-3



"""
    boost1d(spacings::Vector{Float64},frequencies::Vector{Float64};
        eps::Real=24.,thickness::Real=1e-3,tand::Real=0.)

Return analytical1d boost values for given `spacings` and `frequencies`.
"""
function boost1d(spacings::Vector{Float64},frequencies::Vector{Float64};
        eps::Real=24.,thickness::Real=1e-3,tand::Real=0.)
    
    return abs2.(disk_system(frequencies;
        tand=tand,spacings=[spacings;0],disk_thickness=thickness,
        disk_epsilon=eps,num_disk=length(spacings))[2])
end



"""
    ref1d(spacings::Vector{Float64},frequencies::Vector{Float64};
        eps::Real=24.,thickness::Real=1e-3,tand::Real=0.)

Return (complex) analytical1d reflectivity values for given `spacings` and
`frequencies`.
"""
function ref1d(spacings::Vector{Float64},frequencies::Vector{Float64};
        eps::Real=24.,thickness::Real=1e-3,tand::Real=0.)

    return disk_system(frequencies;
        tand=tand,spacings=[spacings;0],disk_thickness=thickness,
        disk_epsilon=eps,num_disk=length(spacings))[1]
end



"""
    getBoost1d(booster::Booster,frequencies::Array{Float64})

Return analytical1d boost values for given `booster` and
    `frequencies`.
"""
function getBoost1d(booster::Booster,frequencies::Array{Float64})
    return boost1d(pos2dist(booster.pos; disk_thickness=booster.thickness),
        frequencies; eps=booster.epsilon,tand=booster.tand,
        thickness=booster.thickness)
end

"""
    getRef1d(booster::Booster,frequencies::Array{Float64})

Return (complex) analytical1d reflectivity values for given `booster` and
`frequencies`.
"""
function getRef1d(booster::Booster,frequencies::Array{Float64})
    return ref1d(pos2dist(booster.pos; disk_thickness=booster.thickness),
        frequencies; eps=booster.epsilon,tand=booster.tand,
        thickness=booster.thickness)
end


"""
    boost3d(spacings::Vector{Float64},frequencies::Vector{Float64};
        eps::Real=24.,tand::Real=0.,thickness::Real=1e-3,R::Real=0.15,
        M::Int=1,L::Int=0,gridwidth::Real=1.0,dx=0.2)

Return 3d boost values for given `spacings` and `frequencies` using Bessel-mode
calculations.

# Arguments
- spacings::Vector{Float64}: Distances between discs.
- frequencies::Vector{Float64}: Frequencies to calculate boost for.
- eps::Real=24.: Relative dielectric constant of disc material.
- tand::Real=0.: Dieletric loss of disc material.
- thickness::Real=1e-3: Thickness of discs.
- R::Real=0.15: Radius of discs.
- M::Int=1: Maximum M value for mode calculations.
- L::Int=0: Maximum L value for mode calculations.
- gridwidth::Real=1.0: Size of discretization grid.
- dx::Real=0.02: Tile size of discretization grid.
"""
function boost3d(spacings::Vector{Float64},frequencies::Vector{Float64};
        eps::Real=24.,tand::Real=0.,thickness::Real=1e-3,R::Real=0.15,
        M::Int=1,L::Int=0,gridwidth::Real=1.0,dx::Real=0.02)

    ndisk = length(spacings)

    eps = eps+1.0im*atan(tand/eps)

    epsilon = ComplexF64[NaN; [isodd(i) ? 1.0 : eps for i in 1:2*ndisk+1]]
    distance = [0.0; [isodd(i) ? spacings[div(i+1,2)] : thickness for i in 1:2*ndisk]; 0.0]

    coords = SeedCoordinateSystem(; X=-gridwidth/2:dx:gridwidth/2,
        Y=-gridwidth/2:dx:gridwidth/2)

    sbdry = SeedSetupBoundaries(coords,diskno=ndisk,distance=distance,epsilon=epsilon)
    modes = SeedModes(coords,ThreeDim=true,Mmax=M,Lmax=L,diskR=R)

    m_reflect = zeros(M*(2*L+1)); m_reflect[L+1] = 1.0

    @everywhere begin
        sbdry = $(sbdry)
        coords = $(coords)
        modes = $(modes)
        m_reflect = $(m_reflect)
        R = $(R)
    end

    boost_total = @sync @distributed (cat_) for f in frequencies
        boost, _ = transformer(sbdry,coords,modes; reflect=m_reflect, prop=propagator,
            diskR=R,f=f)

        sum(abs2.(boost))
    end

    return boost_total
end



"""
    getBoost3d(booster::Booster,frequencies::Array{Float64},
        (M,L,gridwidth,dx)::Tuple{Int,Int,Real,Real}=(1,0,1,0.02))

Return 3d boost values for given `booster` and `frequencies` using Bessel-mode calculations.
See [`boost3d`](@ref).
"""
function getBoost3d(booster::Booster,frequencies::Array{Float64},
        (M,L,gridwidth,dx)::Tuple{Int,Int,Real,Real}=(1,0,1,0.02))

    return boost3d(pos2dist(booster.pos; disk_thickness=booster.thickness),
        frequencies; eps=booster.epsilon,tand=booster.tand,
        thickness=booster.thickness,R=booster.R,
        M=M,L=L,gridwidth=gridwidth,dx=dx)
end



"""
    findpeak1d(frequency::Real,ndisk::Int;
        eps::Real=24.,tand::Real=0.,thickness::Real=1e-3,granularity::Int=1000,
        deviation::Real=0.1)
    
Return the best found spacing using analytical1d that maximizes the boost value at the
given `frequency` for `n` equidistant discs. Search for `granularity` steps between
`(1-deviation)*λ`, `(1+deviation)*λ`.
"""
function findpeak1d(frequency::Real,ndisk::Int;
        eps::Real=24.,tand::Real=0.,thickness::Real=1e-3,granularity::Int=1000,
        deviation::Real=0.1)

    λ = 299792458.0/frequency
    B = zeros(granularity)
    D = range(1-deviation; stop=1+deviation,length=granularity)*λ/2

    for i in eachindex(D)
        B[i] = boost1d(ones(ndisk)*D[i],[frequency];
            eps=eps,tand=tand,thickness=thickness)[1]
    end

    return D[findmax(B)[2]]
end

"""
    findpeak3d(frequency::Real,ndisk::Int;
        eps::Real=24.,tand::Real=0.,thickness::Real=1e-3,granularity::Int=1000,
        deviation::Real=0.1)
    
Return the best found spacing using transformer 3d boost that maximizes the boost value at
the given `frequency` for `n` equidistant discs. Search for `granularity` steps between
`(1-deviation)*λ`, `(1+deviation)*λ`. See [`boost3d`](@ref) for remaining parameters.
"""
function findpeak3d(frequency::Real,n::Int,
        (R,M,L,gridwidth,dx)::Tuple{Real,Int,Int,Real,Real}=(0.15,1,0,1,0.02);
        eps::Real=24.,tand::Real=0.,thickness::Real=1e-3,
        granularity::Int=1000,deviation::Real=0.1)

    λ = 299792458.0/frequency
    B = zeros(granularity)
    D = range(1-deviation; stop=1+deviation,length=granularity)*λ/2

    for i in eachindex(D)
        B[i] = boost3d(ones(n)*D[i],[frequency];
            eps=eps,tand=tand,thickness=thickness,R=R,
            M=M,L=L,gridwidth=gridwidth,dx=dx)[1]
    end

    return D[findmax(B)[2]]
end



"""
    genFreqs(fcenter::Real,fwidth::Real; n::Int=100)

Return `n` equally spaced frequencies from `fcenter-fwidth/2` to
`fcenter+fwidth/2`.
"""
function genFreqs(fcenter::Real,fwidth::Real; n::Int=100)
    return Array(range(fcenter-fwidth/2; stop=fcenter+fwidth/2,length=n))
end

"""
    genFreqs(bounds; n::Int=100)

Return `n` equally spaced frequencies from `bounds[1]` to `bounds[2]`.
"""
function genFreqs(bounds; n::Int=100)
    return Array(range(bounds[1]; stop=bounds[2],length=n))
end

"""
    pos2dist(position::Array{Float64}; disk_thickness::Real=1e-3)

Return distances corresponding to `position`.
"""
function pos2dist(position::Array{Float64}; disk_thickness::Real=1e-3)
    pos = [0; position]
    d = (pos[2:end]-pos[1:end-1])
    d[2:end] .-= disk_thickness
    
    return d
end


"""
    dist2pos(distances::Array{Float64}; disc_thickness::Real=1e-3)

Return position corresponding to `distances`.
"""
function dist2pos(distances::Array{Float64}; disc_thickness::Real=1e-3)
    return [sum(distances[1:i])+(i-1)*disc_thickness for i in 1:length(distances)]
end




"""
    pNorm(x,p=2)

Return `(∑ x^p)^(1/p)`.
"""
function pNorm(x,p=2)
    return sum(@. abs(x)^p)^(1/p)
end


Base.copy(x::T) where T = T([getfield(x, k) for k ∈ fieldnames(T)]...)

function shiftdown!(x::Vector)
    @inbounds for i in length(x)-1:-1:1
        x[i+1] = x[i]
    end
end




function getMag(x::Real)
    return Int(3*round(log(10,x)/3,RoundNearest))
end

function magLabel(mag::Int)
    return get(magnitude_labels,mag,"10^$mag ")
end

const magnitude_labels = Dict{Int,String}(
    -9 => "n",
    -6 => "µ",
    -3 => "m",
     0 => "",
     3 => "k",
     6 => "M",
     9 => "G",
     12 => "T"
)

"""
    e(n::Int,idxs::NTuple{N,Int},T::Type=Float64) where N

    e(n::Int,idx::Int,T::Type=Float64)

Return vector of type `T` with all zeros, except for ones at `idxs`.
"""
function e(n::Int,idxs::Array{<:Integer},T::Type=Float64)
    e_ = zeros(T,n)

    e_[idxs] .= one(T)

    return e_
end

function e(n::Int,idx::Int,T::Type=Float64)
    e_ = zeros(T,n); e_[idx] = one(T)

    return e_
end


# ≽^•⩊•^≼              
cat_(args...) = Base.cat(dims=1,args...)