#backbone functions of the optimisers
import Dates: Second, canonicalize



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



###     information tracing

mutable struct Trace
    x::Array{Float64}
    obj::Float64
    g::Vector{Float64}
    h::Matrix{Float64}
    t::Float64
    T::Float64
end

Trace(x,obj,g) = Trace(x,obj,g,zeros(1,1),0.,0.)

Trace() = Trace([0],0,[0],zeros(1,1),0.,0.)

mutable struct NMTrace
    x::Matrix{Float64}
    obj::Vector{Float64}
    x_::Vector{Float64}
    obj_::Float64
    t::Float64
    T::Float64
end

###     old derivative functions

function getGradient(booster::Booster,hist::Array{State},freqs::Array{Float64},
                        Δx::Float64; mode="double")
    gradient = zeros(Float64,booster.ndisk)

    if mode == "double"
        for i in 1:booster.ndisk
            moveCommand(booster,[(i,Δx)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(i,-2Δx)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(i,Δx)])
            gradient[i] = (hist[2].objvalue - hist[1].objvalue)/2Δx
        end
    else
        for i in 1:booster.ndisk
            moveCommand(booster,[(i,Δx)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(i,-Δx)])
            updateHist!(booster,hist,freqs,objFunction)
            gradient[i] = (hist[2].objvalue - hist[1].objvalue)/Δx
        end
    end

    return gradient
end

function getHessian(booster::Booster,hist::Array{State},freqs::Array{Float64},
                    Δx::Float64; mode="forward")
    hessian = zeros(Float64,booster.ndisk,booster.ndisk)

    if mode == "forward"
        for i in 1:booster.ndisk, j in 1:booster.ndisk
            moveCommand(booster,[(i,Δx),(j,Δx)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(j,-Δx)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(i,-Δx),(j,Δx)])
            updateHist!(booster,hist,freqs,objFunction)
            moveCommand(booster,[(j,-Δx)])
            updateHist!(booster,hist,freqs,objFunction)
            hessian[i,j] = (hist[4].objvalue-hist[3].objvalue-
                            hist[2].objvalue+hist[1].objvalue)/Δx^2
        end
    elseif mode == "backward"
        for i in 1:booster.ndisk, j in 1:booster.ndisk
            moveCommand(booster,[i,0])
        end
    end

    return hessian
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

function printIter(booster::Booster,hist,i::Int,k::Int)
    println("Iter: ",i,", timestamp: ",round(booster.timestamp,digits=3))
    println("Iter finished. Steps: ",k,", Objective value: ",
            round(hist[1].objvalue; digits=3),"\n")
    k == 0 && println("Stuck. Trying to unstuck.\n")
end

function printNMIter(booster::Booster,f::Vector{Float64},i::Int)
    println("Iter: ",i,", timestamp: ",
                            canonicalize(Second(round(Int,booster.timestamp))))
    println("Iter finished. Objective value: ",
                                        round(minimum(f); digits=3),"\n")
end

function analyse(hist,trace::Vector{Trace},freqsplot;
                freqs=nothing,plotting=true,div=5,scale=1e9,ylim=[-0.05e4,3e4])
    tracex = hcat((x->x.x).(trace)...)
    traced = hcat((x->pos2dist(x.x)).(trace)...)
    tracef = (x->x.obj).(trace)
    traceg = hcat((x->x.g).(trace)...)
    traceh = cat((x->x.h).(trace)...;dims=3)

    l = length(trace)
    n = length(tracex[:,1])

    if plotting
        plt1 = plot(freqsplot/scale,boost1d(pos2dist(tracex[:,1]),freqsplot);
                ylim=ylim,label="init",lc="blue",lw=2)

        if div != 0
            for i in 2:maximum([1,l÷div]):(l-1)
                plot!(freqsplot/scale,boost1d(pos2dist(tracex[:,i]),freqsplot);
                        ylim=ylim,label="it.: "*string(i))
            end
        end

        plot!(freqsplot/scale,boost1d(pos2dist(tracex[:,l]),freqsplot);
                ylim=ylim,label="final",lc="red",lw=2)

        if freqs != nothing
            vline!([minimum(freqs),maximum(freqs)]/scale,c="black",linestyle=:dash,
                    label="")
        end
        title!("Boostfactor")
        xlabel!("Frequency [GHz]")
        ylabel!("β²")
        annotate!([(minimum(freqsplot)/scale,0.9*ylim[2],
                    "Final value:\n"*string(round(tracef[l],digits=1)),:left)])

        plt2 = plot(1:l,tracef;legend=false)
        title!("Objective trace")
        xlabel!("Iteration")
        ylabel!("Objective value")

        plt3 = plot(1:l,traced';legend=false)
        title!("Distance trace")
        xlabel!("Iteration")
        ylabel!("d_i")

        plt4 = scatter(1:n,traced[:,l];legend=false)
        title!("Final distances")
        xlabel!("Disk")
        ylabel!("d_i")

        plt5 = plot((x->x.objvalue).(hist[(x->x.objvalue).(hist) .!= 0.])[end:-1:1][1:end];
                    legend=false)
        title!("History")
        xlabel!("Step")
        ylabel!("Objective value")

        display(plt1)
        display(plt2)
        display(plt3)
        display(plt4)
        display(plt5)
    end

    if !plotting
        return tracex, traced, tracef, traceg, traceh
    else
        return plt1, plt2, plt3, plt4, plt5
    end
end

function analyse(hist,trace::Vector{NMTrace},freqsplot;
                        freqs=nothing,
                        plotting=true,
                        div=5,
                        scale=1e9,
                        ylim=[-0.05e4,3e4])
    tracex = hcat((x->x.x[:,1]).(trace)...)
    tracex_ = hcat((x->x.x_).(trace)...)
    traced = hcat((x->pos2dist(x.x[:,1])).(trace)...)
    traced_ = hcat((x->pos2dist(x.x_)).(trace)...)
    tracef = (x->x.obj[1]).(trace)
    tracef_ = (x->x.obj_).(trace)

    l = length(trace)
    n = length(tracex[:,1])

    lh = length(hist[(x->x.objvalue).(hist) .!= 0.])

    histx = hcat((x->x.pos).(hist[lh:-1:1])...)
    histf = (x->x.objvalue).(hist[lh:-1:1])
    histd = hcat((x->pos2dist(x.pos)).(hist[lh:-1:1])...)

    if plotting
        plt1 = plot(freqsplot/scale,boost1d(pos2dist(tracex[:,1]),freqsplot);
                ylim=ylim,label="init",lc="blue",lw=2)

        if div != 0
            for i in 2:maximum([1,l÷div]):(l-1)
                plot!(freqsplot/scale,boost1d(pos2dist(tracex[:,i]),freqsplot);
                        ylim=ylim,label="it.: "*string(i))
            end
        end

        plot!(freqsplot/scale,boost1d(pos2dist(tracex[:,l]),freqsplot);
                ylim=ylim,label="final",lc="red",lw=2)

        if freqs != nothing
            vline!([minimum(freqs),maximum(freqs)]/scale,c="black",linestyle=:dash,
                    label="")
        end
        title!("Boostfactor")
        xlabel!("Frequency [GHz]")
        ylabel!("β²")
        annotate!([(minimum(freqsplot)/scale,0.9*ylim[2],
                    "Final value:\n"*string(round(tracef[l],digits=1)),:left)])

        plt2 = plot(1:l,tracef; legend=false)
        title!("Objective trace best vertex")
        xlabel!("Iteration")
        ylabel!("Objective value")

        plt3 = plot(1:l,traced'; legend=false)
        title!("Distance trace best vertex")
        xlabel!("Iteration")
        ylabel!("d_i")

        plt4 = scatter(1:n,traced[:,l]; legend=false)
        title!("Final distances")
        xlabel!("Disk")
        ylabel!("d_i")

        plt5 = plot(1:lh,histf[1:lh]; legend=false)
        title!("History objective value")
        xlabel!("Step")
        ylabel!("Objective value")

        plt6 = plot(1:lh,histd[:,1:lh]'; legend=false)
        title!("History distances")
        xlabel!("Step")
        ylabel!("d_i")

        display(plt1)
        display(plt2)
        display(plt3)
        display(plt4)
        display(plt5)
        display(plt6)
    end

    if !plotting
        return tracex, tracex_, traced, traced_, tracef, tracef_,
                                                            histx, histf, histd
    else
        return plt1, plt2, plt3, plt4, plt5, plt6
    end
end
