
###     physical state of the booster

export Booster, DevicesType, BoundariesType, AnalyticalBooster, PhysicalBooster, State, Callback, F

unow() = now(UTC)

abstract type Booster end

abstract type DevicesType end
abstract type BoundariesType end

abstract type Trace end

mutable struct AnalyticalBooster <: Booster
    pos::Array{<:Real}
    ndisk::Int
    thickness::Real
    epsilon::Float64
    tand::Number
    vmotor::Real
    maxlength::Real
    timestamp::DateTime
    codetimestamp::DateTime
    summeddistance::Float64

    function AnalyticalBooster(initdist; ndisk=20,τ=1e-3,ϵ=24,tand=0,vmotor=0.1e-3,maxlength=2)
        new(dist2pos(initdist*ones(ndisk)),ndisk,τ,ϵ,tand,vmotor,maxlength,DateTime(0),unow(),0.)
    end

    function AnalyticalBooster(pos,ndisk,thickness,epsilon,tand,vmotor,maxlength,timestamp,codetimestamp,summeddistance)
        new(pos,ndisk,thickness,epsilon,tand,vmotor,maxlength,timestamp,codetimestamp,summeddistance)
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
    summeddistance::Float64

    function PhysicalBooster(devices,initdist; ndisk=20,τ=1e-3,ϵ=24,maxlength=2)
        new(devices,dist2pos(initdist*ones(ndisk)),ndisk,τ,ϵ,maxlength,
            unow(),unow(),0.)
    end

    function PhysicalBooster(
            devices,pos,ndisk,thickness,epsilon,maxlength,timestamp,
            startingtime,summeddistance)

        new(devices,pos,ndisk,thickness,epsilon,maxlength,timestamp,
            startingtime,summeddistance)
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
        new(copy(booster.pos),objvalue,booster.timestamp)        
    end

    function State(booster::PhysicalBooster,objvalue)
        new(copy(booster.pos),objvalue,unow())    # fix???     
    end

    function State(booster,objvalue,timestamp)
        new(copy(booster.pos),objvalue,timestamp)
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