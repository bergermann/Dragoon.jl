###     backbone functions for all optimizers

export getState, initHist, updateHist!

function getState(booster::Booster,
                    freqs::Array{Float64},
                    objFunction::Callback)

    return State(booster,objFunction.func(booster,freqs,objFunction.args))
end

function initHist(booster::Booster,
                    length::Int64,
                    freqs::Array{Float64},
                    objFunction::Callback)

    s0 = State(booster)
    hist = fill(s0,length)
    updateHist!(booster,hist,freqs,objFunction; showtrace=true,force=true)

    return hist
end

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
end



###     control functions for the booster

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

    # booster.summeddistance += sum(abs.(booster.pos-newpos))

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

function move(booster::AnalyticalBooster,newpos::Array{Float64};
        Δt=0,returntrace=false,tracestep=1e-3,additive=false)

    booster.summeddistance += sum(abs.(booster.pos-newpos))
    
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

function printTimes(booster::Booster)
    if hasproperty(booster,:startingtime)
        println("Elapsed movement time:  ",canonicalize(round(booster.timestamp-booster.startingtime,Second)))
    else
        println("Elapsed movement time:  ",canonicalize(round(booster.timestamp-DateTime(0),Second)))
    end

    println("Summed distance:   ",round(booster.summeddistance; digits=3))

    if hasproperty(booster,:codetimestamp)
        println("Elapsed computing time: ",canonicalize(booster.codetimestamp-DateTime(0)))
    end
end

function printTermination(booster::Booster,hist,i::Int,maxiter::Int)
    if i == maxiter
        println("Terminated. Max iterations reached.")
    else
        println("Terminated. ",i," Iterations.")
    end

    println("Final objective value: ",round(hist[1].objvalue; digits=3))
    printTimes(booster)
end
