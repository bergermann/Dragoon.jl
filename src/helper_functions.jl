###     backbone functions for all optimizers

export getState, initHist, updateHist!, move

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

Shift all data one index down, write current state to first entry. Last entry of
`hist` is always lost.
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

    returnß
end



###     control functions for the booster

"""
    move(booster::AnalyticalBooster,newpos::Vector{Tuple{Int64,Float64}};
        Δt=0,returntrace=false,tracestep=1e-3,additive=true)

Set position of analytical `booster` to \"move\" it and update time values.
Iteratively set positions at `index` for entries `(index,value)` of `newpos`.
If `additive`, add `value` to recent position, else overwrite it.
Use `Δt` for additional movement overhead.

[Tracing NYI]
"""
function move(booster::AnalyticalBooster,newpos::Vector{Tuple{Int64,Float64}};
        Δt=0,returntrace=false,tracestep=1e-3,additive=true)
    
    T = zeros(length(newpos))

    for n in newpos
        if additive
            booster.summeddistance += abs(n[2])
        else
            booster.summeddistance += abs(booster.pos[n[1]]-n[2])
        end
    end

    if additive
        for i in 1:length(newpos)
            T[i] = abs(newpos[i][2])/booster.vmotor
            booster.pos[newpos[i][1]] += newpos[i][2]
        end
    else
        for i in 1:length(newpos)
            T[i] = abs(p[i][2]-booster.pos[newpos[i][1]])/booster.vmotor
            booster.pos[newpos[i][1]] = newpos[i][2]
        end
    end

    booster.timestamp += (Δt + maximum(T)) *ₜ Second

    if returntrace
        return 0
    end
end

"""
    move(booster::AnalyticalBooster,newpos::Array{Float64};
        Δt=0,returntrace=false,tracestep=1e-3,additive=false)

Set position of analytical `booster` to \"move\" it and update time values.
If `additive`, position is set to `booster.pos + newpos`, else to `newpos`.
Use `Δt` for additional movement overhead.

[Tracing NYI]
"""
function move(booster::AnalyticalBooster,newpos::Array{Float64};
        Δt=0,returntrace=false,tracestep=1e-3,additive=false)

    if additive
        booster.summeddistance += sum(abs.(newpos))
    else
        booster.summeddistance += sum(abs.(booster.pos-newpos))
    end
    
    if additive
        T1 = maximum(abs.(newpos))/booster.vmotor
        T2 = sum(abs.(newpos))/booster.vmotor

        if returntrace
            trace = zeros(length(booster.pos),T/tracestep)
        end

        booster.pos += newpos
    else
        T1 = maximum(abs.(booster.pos-newpos))/booster.vmotor
        T2 = sum(abs.(booster.pos-newpos))/booster.vmotor

        if returntrace
            trace = zeros(length(booster.pos),T/tracestep)
        end

        booster.pos = copy(newpos)
    end

    booster.timestamp += (Δt + T1) *ₜ Second

    if returntrace
        return trace
    end
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



