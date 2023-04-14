
###     physical state of the booster

abstract type Booster end

mutable struct AnalyticalBooster <: Booster
    pos::Array{<:Real}
    ndisk::Int
    thickness::Float64
    epsilon::Real
    vmotor::Real
    maxlength::Real
    timestamp::Float64
    summedtraveltime::Float64
    codetimestamp

    function AnalyticalBooster(initdist; ndisk=20,τ=1e-3,ϵ=24,vmotor=0.1e-3,maxlength=2)
        new(dist2pos(initdist*ones(ndisk)),ndisk,τ,ϵ,vmotor,maxlength,0.,0.,0.)
    end
end

abstract type DevicesType end

mutable struct PhysicalBooster <: Booster
    devices<:DevicesType
    pos::Array{<:Real}
    ndisk::Int
    thickness::Float64
    epsilon::Real
    maxlength::Real
    summedtraveltime::Float64

    function PhysicalBooster(initdist; ndisk=20,τ=1e-3,ϵ=24,maxlength=2)
        new(dist2pos(initdist*ones(ndisk)),ndisk,τ,ϵ,maxlength,0.)
    end
end

mutable struct State
    pos::Array{Float64}
    objvalue::Float64
    timestamp::Float64

    function State(booster)
        new(booster.pos,0.0,0.0)
    end

    function State(booster,objvalue,timestamp)
        new(booster.pos,objvalue,timestamp)
    end
end
