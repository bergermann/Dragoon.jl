export Booster, State, initHist

#backbone functions

###     physical state of the booster

mutable struct Booster
    pos::Array{Float64}
    ndisk::UInt8
    thickness::Float64
    epsilon::Float64
    vmotor::Float64
    maxlength::Float64
    timestamp::Float64
    summedtraveltime::Float64
    codetimestamp
end

Booster() = Booster(dist2pos(init),20,1e-3,24.,0.1e-3,2.,0.,0.,0.)

mutable struct State
    pos::Array{Float64}
    objvalue::Float64
    timestamp::Float64
end

State(booster::Booster) = State(zeros(Float64,booster.ndisk),0.0,0.0)

function getState(booster::Booster,
                    freqs::Array{Float64},
                    objFunction::Tuple{Function,Vector})

    return State(booster.pos,objFunction[1](booster,freqs,objFunction[2]...),
                                                            booster.timestamp)
end

function initHist(booster::Booster,
                    length::Int64,
                    freqs::Array{Float64},
                    objFunction::Tuple{Function,Vector})

    hist = fill(State(booster),length)
    updateHist!(booster,hist,freqs,objFunction; showtrace=true,force=true)

    return hist
end

function updateHist!(booster::Booster,
                        hist,freqs::Array{Float64},
                        objFunction::Tuple{Function,Vector};
                        showtrace=false,force=false)

    if hist[1].timestamp != booster.timestamp || force
        hist[2:end] = hist[1:end-1]
        hist[1] = getState(booster,freqs,objFunction)
    end

    if showtrace
        println("Objective Value: ",round(hist[1].objvalue; digits=1),
                    ", Timestamp: ",hist[1].timestamp)
    end
end



###     control functions for the booster

function moveCommand(booster::Booster,newpos::Vector{Tuple{Int64,Float64}};
        Δt=0,returntrace=false,tracestep=1e-3,additive=true)
    T = zeros(length(newpos))

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

    booster.timestamp += Δt + maximum(T)
    booster.summedtraveltime += sum(T)

    if returntrace
        return 0
    end
end

function moveCommand(booster::Booster,newpos::Array{Float64};
        Δt=0,returntrace=false,tracestep=1e-3,additive=false)
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

    booster.timestamp += Δt + T1
    booster.summedtraveltime += T2

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
    println("Elapsed movement time:  ",
                canonicalize(Second(round(Int,booster.timestamp))))
    println("Summed movement time:   ",
                canonicalize(Second(round(Int,booster.summedtraveltime))))
    println("Elapsed computing time: ",canonicalize(booster.codetimestamp))
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
