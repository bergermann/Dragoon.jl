###     convenient convenience functions for convenience

export init, boost1d, ref1d, findpeak, genFreqs, pos2dist, dist2pos, pNorm,
        copy, getBoost1d, getRef1d, groupdelay

init = [1.00334, 6.94754, 7.1766, 7.22788, 7.19717,
        7.23776, 7.07746, 7.57173, 7.08019, 7.24657,
        7.21708, 7.18317, 7.13025, 7.2198, 7.45585,
        7.39873, 7.15403, 7.14252, 6.83105, 7.42282]*1e-3
  
        


boost1d(spacs::Vector{Float64},f::Vector{Float64};
        eps::Real=24.,thickness::Real=1e-3,tand::Real=0.) = 
    abs2.(disk_system(f;
        tand=tand,spacings=[spacs;0],disk_thickness=thickness,disk_epsilon=eps,
        num_disk=length(spacs))[2])

ref1d(spacs::Vector{Float64},f::Vector{Float64};
        eps::Real=24.,thickness::Real=1e-3,tand::Real=0.) = 
    disk_system(f; spacings=[spacs;0],disk_thickness=thickness,disk_epsilon=eps,
        tand=tand,num_disk=length(spacs))[1]

function findpeak(f0,n; eps=24.,tand=0.,thickness=1e-3,gran=1000,dev=0.1)
    λ = 299792458.0/f0
    B = zeros(gran)
    D = range(1-dev; stop=1+dev,length=gran)*λ/2

    for i in eachindex(D)
        B[i] = boost1d(ones(n)*D[i],[f0];
            eps=eps,tand=tand,thickness=thickness)[1]
    end

    return D[findmax(B)[2]]
end

function pos2dist(pos::Array{Float64}; thickness=1e-3)
    pos = [0; pos]
    d = (pos[2:end]-pos[1:end-1])
    d[2:end] .-= thickness
    
    return d
end

function dist2pos(dist::Array{Float64}; thickness=1e-3)
    return [sum(dist[1:i])+(i-1)*thickness for i in 1:length(dist)]
end

function pNorm(x; p=2)
    return sum(@. abs(x)^p)^(1/p)
end


Base.copy(x::T) where T = T([getfield(x, k) for k ∈ fieldnames(T)]...)

# function getBoost1d(booster::Booster,freqs::Array{Float64})
#     return boost1d(pos2dist(booster.pos; thickness=booster.thickness),freqs;
#         eps=booster.epsilon,tand=booster.tand,thickness=booster.thickness)
# end

# function getRef1d(booster::Booster,freqs::Array{Float64})
#     return ref1d(pos2dist(booster.pos; thickness=booster.thickness),freqs;
#         eps=booster.epsilon,tand=booster.tand,thickness=booster.thickness)
# end



function shiftdown!(x::Vector)
    @inbounds for i in length(x)-1:-1:1
        x[i+1] = x[i]
    end
end

function groupdelay(x::Vector{ComplexF64},freqs::Vector{Float64})
    ϕ = angle.(x)
    dϕ = ϕ[2:end]-ϕ[1:end-1]

    gd = dϕ ./ (freqs[2:end]-freqs[1:end-1])

    return gd
end