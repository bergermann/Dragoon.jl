###     information tracing and output linesearch

export analyse

import Plots: plot, plot!, scatter, vline!, title!, xlabel!, ylabel!,
        annotate!, display


###     information tracing

mutable struct LSTrace <: Trace
    x::Array{Float64}
    obj::Float64
    g::Vector{Float64}
    h::Matrix{Float64}
    t::DateTime
    T::Float64
end

LSTrace(x,obj,g) = LSTrace(x,obj,g,zeros(1,1),0.,0.)

LSTrace() = LSTrace([0],0,[0],zeros(1,1),0.,0.)



###     output

function printIter(booster::Booster,hist,i::Int,k::Int)
    if hasproperty(booster,:startingtime)
        println("Iter: ",i,", timestamp: ",canonicalize(
            floor(booster.timestamp-booster.startingtime,Second)))
    else
        println("Iter: ",i,", timestamp: ",canonicalize(
            floor(booster.timestamp-DateTime(0),Second)))
    end
    
    println("Iter finished. Steps: ",k,", Objective value: ",
            round(hist[1].objvalue; digits=3),"\n")
            
    k == 0 && println("Stuck. Trying to unstuck.\n")
end

function analyse(hist,trace::Vector{LSTrace},freqsplot;
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

        if freqs !== nothing
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
