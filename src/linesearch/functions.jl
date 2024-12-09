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




"""
    analyse(hist,trace::Vector{LSTrace},freqsplot;
        freqs=nothing,plotting=true,div=5,ylim=[-0.05e4,3e4])

Analyse linesearch trace and create plot output using Analytical1d.

# Arguments
- `hist`: Vector containing history of states.
- `trace::Vector{LSTrace}`:  Trace from linesearch algorithm.
- `freqsplot`: Frequencies to plot on. Don't have to be same as optimized freqs.
- `freqs=nothing`: Optimized frequencies, used to show their bounds.
- `plotting=true`: If false, return subtraces as plain vectors instead.
- `div=5`: Amount of intermediate steps.
- `ylim=[-0.05e4,3e4]`: Manual limit of y-axis.
"""
function analyse(booster,hist,trace::Vector{LSTrace},freqsplot;
        freqs=nothing,plotting=true,div=5,ylim=[-0.05e4,3e4])
    
    tracex = hcat((x->x.x).(trace)...)
    traced = hcat((x->pos2dist(x.x)).(trace)...)
    tracef = (x->x.obj).(trace)
    traceg = hcat((x->x.g).(trace)...)
    traceh = cat((x->x.h).(trace)...;dims=3)

    l = length(trace)
    n = length(tracex[:,1])

    lh = length(hist[(x->x.objvalue).(hist) .!= 0.])

    histx = hcat((x->x.pos).(hist[lh:-1:1])...)
    histf = (x->x.objvalue).(hist[lh:-1:1])
    histd = hcat((x->pos2dist(x.pos)).(hist[lh:-1:1])...)

    mag = getMag(maximum(freqsplot)); scale = 10^mag

    if plotting
        plt1 = plot(freqsplot/scale,
            boost1d(Pos,tracex[:,1],freqsplot; eps=booster.epsilon,tand=booster.tand,
            thickness=booster.thickness);
            ylim=ylim,label="init",lc="blue",lw=2)

        if div != 0
            for i in 2:maximum([1,l÷div]):(l-1)
                plot!(plt1,freqsplot/scale,
                    boost1d(Pos,tracex[:,i],freqsplot;eps=booster.epsilon,
                    tand=booster.tand,thickness=booster.thickness);
                    ylim=ylim,label="it.: "*string(i))
            end
        end

        plot!(plt1,freqsplot/scale,
            boost1d(Pos,tracex[:,l],freqsplot;eps=booster.epsilon,tand=booster.tand,
            thickness=booster.thickness);
            ylim=ylim,label="final",lc="red",lw=2)

        if freqs !== nothing
            vline!(plt1,[minimum(freqs),maximum(freqs)]/scale,c="black",
                linestyle=:dash,label="")
        end

        title!("Boostfactor")
        xlabel!("Frequency [$(magLabel(mag))Hz]")
        ylabel!(L"Power Boost Factor $β^2$")
        annotate!([(minimum(freqsplot)/scale,0.9*ylim[2],
                    "Final value:\n"*string(round(tracef[l],digits=1)),:left)])

        plt2 = plot(1:l,tracef; legend=false)
        title!("Objective trace best vertex")
        xlabel!("Iteration")
        ylabel!(L"Objective value $f$")

        plt3 = plot(1:l,traced'/1e-3; legend=false)
        title!("Distance trace best vertex")
        xlabel!("Iteration")
        ylabel!(L"Distances $d_i$ [mm]")

        plt4 = scatter(1:n,traced[:,l]/1e-3; legend=false)
        title!("Final distances")
        xlabel!("Disk index")
        ylabel!(L"Distances $d_i$ [mm]")

        plt5 = plot(-lh:-1,histf[1:lh]; legend=false)
        title!("History objective value")
        xlabel!("Step index")
        ylabel!(L"Objective value $f$")

        plt6 = plot(-lh:-1,histd[:,1:lh]'/1e-3; legend=false)
        title!("History distances")
        xlabel!("Step index")
        ylabel!(L"Distances $d_i$ [mm]")

        display(plt1)
        display(plt2)
        display(plt3)
        display(plt4)
        display(plt5)
        display(plt6)

        return plt1, plt2, plt3, plt4, plt5, plt6
    else
        return tracex, traced, tracef, traceg, traceh, histx, histf, histd
    end
end
