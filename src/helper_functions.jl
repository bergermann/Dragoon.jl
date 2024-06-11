###     backbone functions for all optimizers

export getState, initHist, updateHist!, move, modb

"""
    getState(booster::Booster,
        freqs::Array{Float64},
        objFunction::Callback)

Pack `booster` information and objective value into State type. See
[`State`](@ref).
"""
function getState(booster::Booster,freqs::Array{Float64},objFunction::Callback)
    return State(booster,objFunction.func(booster,freqs,objFunction.args))
end

"""
    initHist(booster::Booster,
        length::Int64,
        freqs::Array{Float64},
        objFunction::Callback;
        showtrace::Bool=false)

Initialize history vector filled with zeros, perform single update. See
[`updateHist!`](@ref)
"""
function initHist(booster::Booster,
        length::Int64,
        freqs::Array{Float64},
        objFunction::Callback;
        showtrace::Bool=false)

    s0 = State(booster)
    hist = fill(s0,length)
    updateHist!(booster,hist,freqs,objFunction; showtrace=showtrace,force=true)

    return hist
end

"""
    updateHist!(booster::Booster,
        hist::Vector{State},
        freqs::Array{Float64},
        objFunction::Callback;
        showtrace=false,force=false)

Shift all history data one index down, write current state to first index. Last entry of
`hist` is always lost. Update only if current position is not last position, unless
`forced`. Return objective value.
"""
function updateHist!(booster::Booster,
        hist::Vector{State},
        freqs::Array{Float64},
        objFunction::Callback;
        showtrace=false,force=false)

    if hist[1].pos != booster.pos || force
        shiftdown!(hist)
        hist[1] = getState(booster,freqs,objFunction)
    end

    if showtrace
        println("Objective Value: ",round(hist[1].objvalue; digits=1),
                    ", Timestamp: ",hist[1].timestamp)
    end

    return hist[1].objvalue
end



###     control functions for the booster

"""
    move(booster::AnalyticalBooster,newpos::Vector{Tuple{Int64,Float64}};
        Δt=0,tracestep=1e-3,additive=true)

Set position of analytical `booster` to \"move\" it and update time values.
Iteratively set positions at `index` for entries `(index,value)` of `newpos`.
If `additive`, add `value` to recent position, else overwrite it.
Use `Δt` for additional movement overhead. Return trace of movement [NYI].
"""
function move(booster::AnalyticalBooster,newpos::Vector{Tuple{Int64,Float64}};
        Δt=0,tracestep=1e-3,additive=true)
    
    newpos_ = copy(booster.pos)

    if additive
        for i in eachindex(newpos)
            newpos_[newpos[i][1]] += newpos[i][2]
        end
    else
        for i in eachindex(newpos)
            newpos_[newpos[i][1]] = newpos[i][2]
        end
    end

    return move(booster,newpos_; Δt=Δt,tracestep=tracestep,additive=false)
end

"""
    move(booster::AnalyticalBooster,newpos::Array{Float64};
        Δt=0,tracestep=1e-3,additive=false)

Set position of analytical `booster` to \"move\" it and update time values.
If `additive`, position is set to `booster.pos + newpos`, else to `newpos`.
Use `Δt` for additional movement overhead. Return trace of movement [NYI].
"""
function move(booster::AnalyticalBooster,newpos::Vector{Float64};
        Δt=0,tracestep=1e-3,additive=false)
    
    if additive
        newpos .+= booster.pos
    end

    trace = zeros(length(booster.pos),ceil(Int,maximum(abs.(booster.pos-newpos))/tracestep))

    d = pos2dist(newpos)
    if booster.wavelength != 0
        # d[d .<= 0] .+= booster.wavelength
        @. d = modp(d,booster.wavelength,2*booster.wavelength)
    else
        d .= max.(d,0)
    end
    newpos .= dist2pos(d)

    booster.summeddistance += sum(abs.(booster.pos-newpos))
    booster.timestamp += (Δt + maximum(abs.(booster.pos-newpos))/booster.vmotor) *ₜ Second

    booster.pos = copy(newpos)

    return trace
end

function modb(booster::Booster,newpos::Vector{Float64})
    d = pos2dist(newpos)

    if booster.wavelength != 0
        # d[d .<= 0] .+= booster.wavelength
        @. d = modp(d,booster.wavelength,2*booster.wavelength)
    end

    d[d .<= booster.mindist] .+= booster.wavelength

    return dist2pos(d)
end



###     time calculations

function travelTime(pos1,pos2; speed=0.01)
    return maximum(abs.(pos1-pos2))/speed
end

function totalTravelTime(pos1,pos2; speed=0.01)
    return sum(abs.(pos1-pos2))/speed
end



###     information output functions

function printTimes(booster::Booster,showtrace::Bool)
    ttotal = 0.

    if hasproperty(booster,:startingtime)
        ttotal = round(booster.timestamp-
            booster.startingtime,Second)
        showtrace && println("Elapsed movement time:  ",
            canonicalize(ttotal))
    else
        ttotal = round(booster.timestamp-DateTime(0),Second)
        showtrace && println("Elapsed movement time:  ",
            canonicalize(ttotal))
    end

    sumdist = round(booster.summeddistance; digits=3)
    showtrace && println("Summed distance:   ",sumdist)

    tcomp = 0.
    if hasproperty(booster,:codetimestamp)
        tcomp = booster.codetimestamp-DateTime(0)
        showtrace && println("Elapsed computing time: ",
            canonicalize(tcomp))
    end

    return ttotal, sumdist, tcomp
end

function printTermination(booster::Booster,hist,i::Int,maxiter::Int,showtrace::Bool)
    if i == maxiter
        showtrace && println("Terminated. Max iterations reached.")
    else
        showtrace && println("Terminated. ",i," Iterations.")
    end

    showtrace && println("Final objective value: ",round(hist[1].objvalue; digits=3))
    
    return hist[1].objvalue, printTimes(booster,showtrace)
end




function setTimes!(booster::Booster,reset::Bool)
    if hasproperty(booster, :startingtime) && reset
        t0 = DateTime(0)

        showtrace && println("Resetting starting time.")

        booster.startingtime = unow()
        booster.timestamp = unow()
    elseif hasproperty(booster, :codetimestamp)
        t0 = unow()

        if reset
            booster.timestamp = DateTime(0)
        end
    end

    return t0
end

function updateTimeStamp!(booster::Booster,name::Symbol,reset::Bool,t0::DateTime)
    if hasproperty(booster, name)
        if reset
            setfield!(booster,name,DateTime(0))
        end

        setfield!(booster,name,getfield(booster,name)+(unow()-t0))
    end

    return
end



