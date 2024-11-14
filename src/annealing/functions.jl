




function printSAIter(booster::Booster,obj::Float64,objsol::Float64,τ::Float64,iter::Int)
    if hasproperty(booster,:startingtime)
        println("Iter: ",iter,", timestamp: ",canonicalize(
            round(booster.timestamp-booster.startingtime,Second)))
    else
        println("Iter: ",iter,", timestamp: ",canonicalize(
            round(booster.timestamp-DateTime(0),Second)))
    end

    println("Iter finished. Objective value current:  ",round(obj; digits=3))
    println("               Objective value solution: ",round(objsol; digits=3))
    println("               Temperature:              ",round(τ; digits=3),"\n")
end



mutable struct SATrace <: Trace
    x::Array{Float64}
    obj::Float64
    xsol::Array{Float64}
    objsol::Float64
    τ::Float64
    iter::Int
    t::DateTime
    T::Float64

    function SATrace(x,obj,xsol,objsol,τ,iter,t,T)
        new(x,obj,xsol,objsol,τ,iter,t,T)
    end

    function SATrace()
        new([0],0,0,0,DateTime(0),0)
    end

    function SATrace(x,obj,τ,iter)
        new(x,obj,τ,iter,DateTime(0),0)
    end    
end



function findNeighbour(booster::Booster,rmax::Float64)
    x_ = 2*rand(Float64,booster.ndisk) .- 1
    x_ /= pNorm(x_)

    return rmax*rand()*x_
end

function thermal(objx::Float64,objy::Float64,T::Float64)
    return exp((objx-objy)/T)
end



"""
    analyse(hist,trace::Vector{SATrace},freqsplot;
        freqs=nothing,plotting=true,div=5,ylim=[-0.05e4,3e4])

Analyse simulated annealing trace and create plot output using Analytical1d.

# Arguments
- `hist`: Vector containing history of states.
- `trace::Vector{LSTrace}`:  Trace from linesearch algorithm.
- `freqsplot`: Frequencies to plot on. Don't have to be same as optimized freqs.
- `freqs=nothing`: Optimized frequencies, used to show their bounds.
- `plotting=true`: If false, return subtraces as plain vectors instead.
- `div=5`: Amount of intermediate steps.
- `ylim=[-0.05e4,3e4]`: Manual limit of y-axis.
"""
function analyse(booster,hist,trace::Vector{SATrace},freqsplot;
        freqs=nothing,plotting=true,div=5,ylim=[-0.05e4,3e4])

    tracex = hcat((x -> x.x).(trace)...)
    traced = hcat((x -> pos2dist(x.x)).(trace)...)

    tracexsol = hcat((x -> x.xsol).(trace)...)
    tracedsol = hcat((x -> pos2dist(x.xsol)).(trace)...)

    traceobj = (x -> x.obj).(trace)
    traceobjsol = (x -> x.objsol).(trace)

    tracet = (x -> x.τ).(trace)

    lh = length(hist[(x->x.objvalue).(hist) .!= 0.])

    histx = hcat((x->x.pos).(hist[lh:-1:1])...)
    histf = (x->x.objvalue).(hist[lh:-1:1])
    histd = hcat((x->pos2dist(x.pos)).(hist[lh:-1:1])...)

    l = length(trace)
    n = length(tracex[:, 1])

    mag = getMag(maximum(freqsplot)); scale = 10^mag

    if plotting
        plt1 = plot(freqsplot/scale,
            boost1d(Pos,tracex[:,1],freqsplot;eps=booster.epsilon,
            tand=booster.tand,thickness=booster.thickness);
            ylim=ylim,label="init",lc="blue",lw=2)
        if div != 0
            for i in 2:maximum([1, l ÷ div]):(l-1)
                plot!(freqsplot/scale,
                    boost1d(Pos,tracex[:,i],freqsplot;eps=booster.epsilon,
                    tand=booster.tand,thickness=booster.thickness);
                    ylim=ylim,label="it.: "*string(i))
            end
        end

        plot!(freqsplot/scale,
            boost1d(Pos,tracex[:,l],freqsplot;eps=booster.epsilon,
            tand=booster.tand,thickness=booster.thickness);
            ylim=ylim,label="final",lc="red",lw=2)

        if freqs !== nothing
            vline!([minimum(freqs),maximum(freqs)]/scale,c="black",linestyle=:dash,
                label="")
        end

        title!("Boostfactor")
        xlabel!("Frequency [$(magLabel(mag))Hz]")
        ylabel!(L"Power Boost Factor $β^2$")
        annotate!([(minimum(freqsplot)/scale,0.9*ylim[2],
                    "Final value:\n"*string(round(traceobjsol[l],digits=1)),:left)])

        plt2 = plot(1:l,traceobjsol; legend=false)
        title!("Objective trace best vertex")
        xlabel!("Iteration")
        ylabel!(L"Objective value $f$")

        plt3 = plot(1:l,tracedsol'; legend=false)
        title!("Distance trace best vertex")
        xlabel!("Iteration")
        ylabel!(L"Distances $d_i$ [mm]")

        plt4 = scatter(1:n,tracedsol[:,l]; legend=false)
        title!("Final distances")
        xlabel!("Disk index")
        ylabel!(L"Distances $d_i$ [mm]")

        plt5 = plot(1:l,traceobj; legend=false)
        title!("Thermal objective trace")
        xlabel!("Iteration")
        ylabel!(L"Objective value $f$")

        plt6 = plot(1:l,traced'; legend=false)
        title!("Thermal distance trace")
        xlabel!("Iteration")
        ylabel!(L"Distances $d_i$ [mm]")

        plt7 = plot(-lh:-1,histf[1:lh]; legend=false)
        title!("History objective value")
        xlabel!("Step index")
        ylabel!(L"Objective value $f$")

        plt8 = plot(-lh:-1,histd[:,1:lh]'; legend=false)
        title!("History distances")
        xlabel!("Step index")
        ylabel!(L"Distances $d_i$ [mm]")

        display(plt1)
        display(plt2)
        display(plt3)
        display(plt4)
        display(plt5)
        display(plt6)
        display(plt7)
        display(plt8)

        return plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8
    else
        return tracex, traced, tracexsol, tracedsol, traceobj, traceobjsol, tracet,
            histx, histf, histd
    end
end