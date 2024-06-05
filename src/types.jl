
###     physical state of the booster

export Booster,
    AnalyticalBooster, PhysicalBooster,
    DevicesType, BoundariesType, 
    State, Callback, F

unow() = now(UTC)

"""
    Booster

Abstract supertype of boosters.
"""
abstract type Booster end

abstract type DevicesType end
abstract type BoundariesType end

abstract type Trace end

"""
    AnalyticalBooster <: Booster

Booster type for analytical calculations.

- `AnalyticalBooster(initdist::Real; ndisk=20,ϵ=24,tand=0,τ=1e-3,R=0.15,
        vmotor=0.1e-3,maxlength=2)`

- `AnalyticalBooster(initpos::Vector; ϵ=24,tand=0,τ=1e-3,R=0.15,
        vmotor=0.1e-3,maxlength=2)`
    
- `AnalyticalBooster(pos,ndisk,epsilon,tand,thickness,R,
        vmotor,maxlength,timestamp,codetimestamp,summeddistance,wavelength)`
"""
mutable struct AnalyticalBooster <: Booster
    "Disc positions."
    pos::Array{<:AbstractFloat}
    "Amount of discs."
    ndisk::Integer
    "Disc relative dielectric constant."
    epsilon::AbstractFloat
    "Disc dielectric loss angle."
    tand::AbstractFloat
    "Disc thickness τ in meter."
    thickness::AbstractFloat
    "Disc radius in meter."
    R::Real
    "Velocity of hypothetical motor."
    vmotor::Real
    "Maximum allowed booster length in meter."
    maxlength::Real
    "Summed movement time of discs."
    timestamp::DateTime
    "Elapsed runtime."
    codetimestamp::DateTime
    "Summed movement distance of discs."
    summeddistance::AbstractFloat
    "Wavelength (of center optimization frequency) for reflection at zero."
    wavelength::AbstractFloat

    function AnalyticalBooster(initdist::Real; ndisk=20,ϵ=24,tand=0,τ=1e-3,R=0.15,
            vmotor=0.1e-3,maxlength=2)
        new(dist2pos(initdist*ones(ndisk)),ndisk,ϵ,tand,τ,R,
            vmotor,maxlength,DateTime(0),unow(),0,0)
    end

    function AnalyticalBooster(initpos::Vector; ϵ=24,tand=0,τ=1e-3,R=0.15,
            vmotor=0.1e-3,maxlength=2)
        new(initpos,length(initpos),ϵ,tand,τ,R,vmotor,maxlength,DateTime(0),DateTime(0),0,0)
    end

    function AnalyticalBooster(pos,ndisk,epsilon,tand,thickness,R,
            vmotor,maxlength,timestamp,codetimestamp,summeddistance,wavelength)
        new(pos,ndisk,epsilon,tand,thickness,R,
            vmotor,maxlength,timestamp,codetimestamp,summeddistance,wavelength)
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


"""
    State

Type for booster states.

- `State(booster)`
- `State(booster::AnalyticalBooster,objvalue)`
- `State(booster::PhysicalBooster,objvalue)`
- `State(booster,objvalue,timestamp)`
"""
mutable struct State
    "Disc positions."
    pos::Array{Float64}
    "Objective value at current position."
    objvalue::Float64
    "Booster timestamp."
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

const States::Type = Vector{State}

"""
    Callback

Type combining a function with flexible, extra arguments for a later call.
NOTE: args are not necessarily the only arguements of func!

Alias: `F`

- `Callback(func,args)`
- `Callback(func)`
"""
mutable struct Callback
    "Function to be called."
    func::Function
    "Arguments to add on function call."
    args::Tuple

    function Callback(func,args)
        new(func,args)
    end

    function Callback(func)
        new(func,())
    end
end

"Alias for `Callback`."
const F = Callback