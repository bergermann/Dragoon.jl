###     convenient convenience functions for convenience

export spacings_stefan,
    boost1d, ref1d, getBoost1d, getRef1d,
    findpeak, genFreqs, 
    pos2dist, dist2pos,
    pNorm, copy,

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
    findpeak(frequency::Real,n::Int;
        eps::Real=24.,tand::Real=0.,thickness::Real=1e-3,granularity::Int=1000,
        deviation::Real=0.1)
    
Return the best found spacing that maximizes the boost value at the given
`frequency` for `n` equidistant discs. Search for `granularity` steps between
`(1-deviation)*λ`, `(1+deviation)*λ`.
"""
function findpeak(frequency::Real,n::Int;
        eps::Real=24.,tand::Real=0.,thickness::Real=1e-3,granularity::Int=1000,
        deviation::Real=0.1)

    λ = 299792458.0/frequency
    B = zeros(gran)
    D = range(1-deviation; stop=1+deviation,length=granularity)*λ/2

    for i in eachindex(D)
        B[i] = boost1d(ones(n)*D[i],[frequency];
            eps=eps,tand=tand,thickness=thickness)[1]
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
    d = (position[2:end]-position[1:end-1])
    d[2:end] .-= disk_thickness
    
    return d
end


"""
    dist2pos(distances::Array{Float64}; disc_thickness::Real=1e-3)

Return position corresponding to `distances`.
"""
function dist2pos(dist::Array{Float64}; disc_thickness::Real=1e-3)
    return [sum(dist[1:i])+(i-1)*disc_thickness for i in 1:length(distances)]
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

magnitude_labels = Dict{Int,String}(
    -9 => "n",
    -6 => "µ",
    -3 => "m",
     0 => "",
     3 => "k",
     6 => "M",
     9 => "G",
     12 => "T"
)