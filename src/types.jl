
###     physical state of the booster

export Booster, DevicesType, BoundariesType, AnalyticalBooster, PhysicalBooster, State, Callback, F

unow() = now(UTC)

abstract type Booster end

abstract type DevicesType end
abstract type BoundariesType end

mutable struct AnalyticalBooster <: Booster
    pos::Array{<:Real}
    ndisk::Int
    thickness::Real
    epsilon::Float64
    vmotor::Real
    maxlength::Real
    timestamp::DateTime
    codetimestamp::DateTime
    summedtraveltime::Float64

    function AnalyticalBooster(initdist; ndisk=20,τ=1e-3,ϵ=24,vmotor=0.1e-3,maxlength=2)
        new(dist2pos(initdist*ones(ndisk)),ndisk,τ,ϵ,vmotor,maxlength,unow(),unow(),0.)
    end

    function AnalyticalBooster(pos,ndisk,thickness,epsilon,vmotor,maxlength,timestamp,codetimestamp,summedtraveltime)
        new(pos,ndisk,thickness,epsilon,vmotor,maxlength,timestamp,codetimestamp,summedtraveltime)
    end
end

mutable struct PhysicalBooster <: Booster
    devices::DevicesType
    pos::Array{<:Real}
    # bounds::Array{BoundariesType}
    ndisk::Int
    thickness::Real
    epsilon::Real
    maxlength::Real
    timestamp::DateTime
    startingtime::DateTime
    summedtraveltime::Float64

    function PhysicalBooster(devices,initdist; ndisk=20,τ=1e-3,ϵ=24,maxlength=2)
        new(devices,dist2pos(initdist*ones(ndisk)),ndisk,τ,ϵ,maxlength,
            unow(),unow(),0.)
    end

    function PhysicalBooster(
            devices,pos,ndisk,thickness,epsilon,maxlength,timestamp,
            startingtime,summedtraveltime)

        new(devices,pos,ndisk,thickness,epsilon,maxlength,timestamp,
            startingtime,summedtraveltime)
    end
end


mutable struct State
    pos::Array{Float64}
    objvalue::Float64
    timestamp::DateTime

    function State(booster)
        new(zeros(Float64,length(booster.pos)),0.0,DateTime(0))
    end

    function State(booster::AnalyticalBooster,objvalue)
        new(booster.pos,objvalue,booster.timestamp)        
    end

    function State(booster::PhysicalBooster,objvalue)
        new(booster.pos,objvalue,unow())    # fix???     
    end

    function State(booster,objvalue,timestamp)
        new(booster.pos,objvalue,timestamp)
    end
end

mutable struct Callback
    func::Function
    args::Tuple

    function Callback(func,args)
        new(func,args)
    end

    function Callback(func)
        new(func,())
    end
end

const F = Callback