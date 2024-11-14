###     information tracing and output linesearch

export analyse

import Plots: plot, plot!, scatter, vline!, title!, xlabel!, ylabel!,
        annotate!, display


###     information tracing

mutable struct ATrace <: Trace
    x::Array{Float64}
    obj::Float64
    g::Vector{Float64}
    t::DateTime
    T::Float64
end

ATrace(x,obj,g) = ATrace(x,obj,g,0.,0.)

ATrace() = ATrace([0],0,[0],0.,0.)



###     output

function printAIter(booster::Booster,hist,i::Int)
    if hasproperty(booster,:startingtime)
        println("Iter: ",i,", timestamp: ",canonicalize(
            floor(booster.timestamp-booster.startingtime,Second)))
    else
        println("Iter: ",i,", timestamp: ",canonicalize(
            floor(booster.timestamp-DateTime(0),Second)))
    end
    
    println("Iter finished. Objective value: ",round(hist[1].objvalue; digits=3),"\n")
            
    # println("Stuck. Trying to unstuck.\n")

    return
end




"""
analyse(booster,hist,trace::Vector{ATrace},freqsplot;
    freqs=nothing,plotting=true,div=5,ylim=[-0.05e4,3e4])

Analyse linesearch trace and create plot output using Analytical1d.

# Arguments
- `hist`: Vector containing history of states.
- `trace::Vector{ATrace}`: Trace from adam algorithm.
- `freqsplot`: Frequencies to plot on. Don't have to be same as optimized freqs.
- `freqs=nothing`: Optimized frequencies, used to show their bounds.
- `plotting=true`: If false, return subtraces as plain vectors instead.
- `div=5`: Amount of intermediate steps.
- `ylim=[-0.05e4,3e4]`: Manual limit of y-axis.
"""
function analyse(booster,hist,trace::Vector{ATrace},freqsplot;
        freqs=nothing,plotting=true,div=5,ylim=[-0.05e4,3e4])
    
    tracex = hcat((x->x.x).(trace)...)
    traced = hcat((x->pos2dist(x.x)).(trace)...)
    tracef = (x->x.obj).(trace)
    traceg = hcat((x->x.g).(trace)...)
    traceh = cat((x->x.h).(trace)...;dims=3)

    l = length(trace)
    n = length(tracex[:,1])

    mag = getMag(maximum(freqsplot)); scale = 10^mag

    if plotting
        plt1 = plot(freqsplot/scale,boost1d(pos2dist(tracex[:,1]),freqsplot);
                ylim=ylim,label="init",lc="blue",lw=2)

        if div != 0
            for i in 2:maximum([1,l÷div]):(l-1)
                plot!(plt1,freqsplot/scale,
                    boost1d(pos2dist(tracex[:,i]),freqsplot);
                    ylim=ylim,label="it.: "*string(i))
            end
        end

        plot!(plt1,freqsplot/scale,boost1d(pos2dist(tracex[:,l]),freqsplot);
            ylim=ylim,label="final",lc="red",lw=2)

        if freqs !== nothing
            vline!(plt1,[minimum(freqs),maximum(freqs)]/scale,c="black",
                linestyle=:dash,label="")
        end

        title!(plt1,"Boostfactor")
        xlabel!(plt1,"Frequency [$(magLabel(mag))Hz]")
        ylabel!(plt1,"β²")
        annotate!(plt1,[(minimum(freqsplot)/scale,0.9*ylim[2],
                    "Final value:\n"*string(round(tracef[l],digits=1)),:left)])

        plt2 = plot(1:l,tracef;legend=false)
        title!(plt2,"Objective trace")
        xlabel!(plt2,"Iteration")
        ylabel!(plt2,"Objective value")

        plt3 = plot(1:l,traced';legend=false)
        title!(plt3,"Distance trace")
        xlabel!(plt3,"Iteration")
        ylabel!(plt3,"d_i")

        plt4 = scatter(1:n,traced[:,l];legend=false)
        title!(plt4,"Final distances")
        xlabel!(plt4,"Disk")
        ylabel!(plt4,"d_i")

        x_obj = (x->x.objvalue).(hist[(x->x.objvalue).(hist) .!= 0.])
        plt5 = plot(x_obj[end-1:-1:1]; legend=false)
        title!(plt5,"History")
        xlabel!(plt5,"Step")
        ylabel!(plt5,"Objective value")

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
