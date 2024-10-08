###     information tracing and output for nelder mead

export analyse

import Plots: plot, plot!, scatter, vline!, title!, xlabel!, ylabel!,
        annotate!, display


mutable struct NMTrace <: Trace
    x::Matrix{Float64}
    obj::Vector{Float64}
    x_::Vector{Float64}
    obj_::Float64
    t::DateTime
    T::Float64
end

function printNMIter(booster::Booster,f::Vector{Float64},i::Int)
    if hasproperty(booster,:startingtime)
        println("Iter: ",i,", timestamp: ",canonicalize(
            round(booster.timestamp-booster.startingtime,Second)))
    else
        println("Iter: ",i,", timestamp: ",canonicalize(
            round(booster.timestamp-DateTime(0),Second)))
    end

    println("Iter finished. Objective value: ",round(minimum(f); digits=3),"\n")
end



"""
    analyse(hist,trace::Vector{NMTrace},freqsplot;
    freqs=nothing,
    plotting=true,
    div=5,
    ylim=[-0.05e4,3e4])

Analyse neldermead trace and create plot output using Analytical1d.

# Arguments
- `hist`: Vector containing history of states.
- `trace::Vector{LSTrace}`:  Trace from linesearch algorithm.
- `freqsplot`: Frequencies to plot on. Don't have to be same as optimized freqs.
- `freqs=nothing`: Optimized frequencies, used to show their bounds.
- `plotting=true`: If false, return subtraces as plain vectors instead.
- `div=5`: Amount of intermediate steps.
- `ylim=[-0.05e4,3e4]`: Manual limit of y-axis.
"""
function analyse(booster,hist,trace::Vector{NMTrace},freqsplot;
        freqs=nothing,plotting=true,div=5,ylim=[-0.05e4,3e4])
    
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
    
    mag = getMag(maximum(freqsplot)); scale = 10^mag

    if plotting
        plt1 = plot(freqsplot/scale,
            boost1d(pos2dist(tracex[:,1]),freqsplot; eps=booster.epsilon,tand=booster.tand,
                thickness=booster.thickness);
                ylim=ylim,label="init",lc="blue",lw=2)

        if div != 0
            for i in 2:maximum([1,l÷div]):(l-1)
                plot!(freqsplot/scale,
                    boost1d(pos2dist(tracex[:,i]),freqsplot;eps=booster.epsilon,
                        tand=booster.tand,thickness=booster.thickness);
                        ylim=ylim,label="it.: "*string(i))
            end
        end

        plot!(freqsplot/scale,
            boost1d(pos2dist(tracex[:,l]),freqsplot;eps=booster.epsilon,tand=booster.tand,
                thickness=booster.thickness);
                ylim=ylim,label="final",lc="red",lw=2)

        if freqs !== nothing
            vline!([minimum(freqs),maximum(freqs)]/scale,c="black",linestyle=:dash,
                    label="")
        end
        title!("Boostfactor")
        xlabel!("Frequency [$(magLabel(mag))Hz]")
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

        return plt1, plt2, plt3, plt4, plt5, plt6
    else
        return tracex, tracex_, traced, traced_, tracef, tracef_,
            histx, histf, histd
    end
end


"""
    getSimplexSize(x::Matrix{Float64},f::Vector{Float64})

Return furthest vertex distance from best vertex.
"""
function getSimplexSize(x::Matrix{Float64},f::Vector{Float64})
    idx = argmin(f)

    s = 0

    for i in axes(x,2)
        if i == idx; continue; end

        d = pNorm(x[:,i]-x[:,idx])

        if d > s
            s = d
        end
    end

    return s
end



"""
    getSimplexInnerSize(x::Matrix{Float64})

Return closest vertex distance from simplex centroid.
"""
function getSimplexInnerSize(x::Matrix{Float64})
    d = Inf64

    x_ = sum(x; dims=2)/size(x,2)

    for i in axes(x,2)
        d_ = pNorm(x[:,i]-x_)

        if d_ < d
            d = d_
        end
    end

    return d
end