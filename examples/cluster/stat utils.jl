using Statistics
using LaTeXStrings

function plotTerm(T,sigx; threshold=-10_000,maxt=Inf,maxdist=Inf,
        annpos=[(0.95,0.2),(0.95,0.2),(0.95,0.2),(0.95,0.2)],ylims=(-15,0.1))
    
    obj = T[:,1,:]
    ttot = T[:,2,:]
    dist = T[:,3,:]

    lx1 = min(maximum(ttot),maxt)

    p1 = scatter(ttot,obj/1e3; legend=false,markersize=2,layout=grid(2,2),
        yflip=true)

    xlims!(p1,(0,lx1*1.01)); ylims!(p1,ylims)

    for s in eachindex(sigx)
        N = size(T,1); nsuccess = sum(T[:,1,s] .<= threshold)
        succ = nsuccess/N; avg = round(Int,mean(T[:,1,s]))

        annotate!(p1[s],lx1*annpos[s][1],ylims[1]*annpos[s][2],
            (
                L"\sigma_x="*"$(round(Int,sigx[s]/1e-6)) µm\n"*
                L"\eta="*"$succ\n"*
                L"\overline{f_{\beta,F}}="*"$avg",
                8,:black,:right
            ))
    end

    plot!(p1[1]; xformatter=:none,bottom_margin=(-3.5,:mm),right_margin=(0.,:mm))
    plot!(p1[2]; xformatter=:none, yformatter=:none,bottom_margin=(-3.5,:mm),
        left_margin=(-3.5,:mm))
    plot!(p1[3]; top_margin=(0.,:mm),xlabel=" ",ylabel=" ")
    plot!(p1[4]; yformatter=:none,left_margin=(-3.5,:mm))

    annotate!(p1[3],lx1,3,(L"\sum\Delta t"*" [s]",12,:black,:center))
    annotate!(p1[3],-0.175*lx1,-16,(L"f_{\beta,F}\times 10^3",12,90.,:center,:black))

    display(p1)

    
    lx2 = min(maximum(dist)/1e-3,maxdist)

    p2 = scatter(dist/1e-3,obj/1e3; legend=false,markersize=2,layout=grid(2,2),
        yflip=true)

    xlims!(p2,(0,lx2*1.01)); ylims!(p2,ylims)

    for s in eachindex(sigx)
        N = size(T,1); nsuccess = sum(T[:,1,s] .<= threshold)
        succ = nsuccess/N; avg = round(Int,mean(T[:,1,s]))

        annotate!(p2[s],lx2*annpos[s][1],ylims[1]*annpos[s][2],
            (
                L"\sigma_x="*"$(round(Int,sigx[s]/1e-6)) µm\n"*
                L"\eta="*"$succ\n"*
                L"\overline{f_{\beta,F}}="*"$avg",
                8,:black,:right
            ))
    end

    plot!(p2[1]; xformatter=:none,bottom_margin=(-3.5,:mm),right_margin=(0.,:mm))
    plot!(p2[2]; xformatter=:none, yformatter=:none,bottom_margin=(-3.5,:mm),
        left_margin=(-3.5,:mm))
    plot!(p2[3]; top_margin=(0.,:mm),xlabel=" ",ylabel=" ")
    plot!(p2[4]; yformatter=:none,left_margin=(-3.5,:mm))

    annotate!(p2[3],lx2,3,(L"\sum\Delta X"*" [mm]",12,:center,:black))
    annotate!(p2[3],-0.175*lx2,-16,(L"f_{\beta,F}\times 10^3",12,90.,:center,:black))

    display(p2)

    return p1, p2
end


import Dragoon: NMTrace, analyse

function analyse(hist, trace::Vector{NMTrace}, freqsplot;
        freqs=nothing,
        plotting=true,
        div=5,
        scale=1e9,
        ylim=[-0.05e4, 3e4],
        freqticks=[21.95,22.0,22.05,22.1])

    tracex = hcat((x -> x.x[:, 1]).(trace)...)
    tracex_ = hcat((x -> x.x_).(trace)...)
    traced = hcat((x -> pos2dist(x.x[:, 1])).(trace)...)
    traced_ = hcat((x -> pos2dist(x.x_)).(trace)...)
    tracef = (x -> x.obj[1]).(trace)/1e3
    tracef_ = (x -> x.obj_).(trace)/1e3

    l = length(trace)
    n = length(tracex[:, 1])

    lh = length(hist[(x->x.objvalue).(hist).!=0.0])

    histx = hcat((x -> x.pos).(hist[lh:-1:1])...)
    histf = (x -> x.objvalue).(hist[lh:-1:1])/1e3
    histd = hcat((x -> pos2dist(x.pos)).(hist[lh:-1:1])...)

    if plotting
        plt1 = plot(freqsplot / scale, boost1d(pos2dist(tracex[:, 1]), freqsplot)/1e3;
            ylim=ylim/1e3, label="init", lc="blue", lw=2)

        if div != 0
            for i in 2:maximum([1, l ÷ div]):(l-1)
                plot!(freqsplot / scale, boost1d(pos2dist(tracex[:, i]), freqsplot)/1e3;
                    ylim=ylim/1e3, label="it.: " * string(i))
            end
        end

        plot!(freqsplot / scale, boost1d(pos2dist(tracex[:, l]), freqsplot)/1e3;
            ylim=ylim/1e3, label="final", lc="red", lw=2)

        if !isnothing(freqs)
            vline!([minimum(freqs), maximum(freqs)] / scale, c="black", linestyle=:dash,
                label="")
        end
        # title!("Boostfactor")
        xticks!(freqticks)
        xlabel!("Frequency [GHz]")
        ylabel!("Power Boost Factor "*L"\beta^2")
        annotate!([(minimum(freqsplot) / scale, 0.9 * ylim[2]/1e3,
            L"f_{\beta,F}:" * string(round(tracef[l]*1e3, digits=0)), :left)])
        annotate!((minimum(freqsplot)+-(extrema(freqsplot)...)*0.1) / scale,
            -0.05*ylim[2]/1e3,(L"\times 10^3",10,:left))

        plt2 = plot(1:l, tracef; legend=false)
        # title!("Objective trace best vertex")
        xlabel!("Iteration")
        ylabel!("Objective Value "*L"f_{\beta}")
        annotate!(-0.1*length(tracef),-(extrema(tracef)...)*0.95,
            (L"\times 10^3",10,:left))

        plt3 = plot(1:l, traced'/1e-3; legend=false)
        # title!("Distance trace best vertex")
        xlabel!("Iteration")
        ylabel!("Distances "*L"d_i"*" [mm]")

        plt4 = scatter(1:n, traced[:, l]/1e-3; legend=false)
        # title!("Final distances")
        xlabel!("Disk")
        ylabel!("Distances "*L"d_i"*" [mm]")

        plt5 = plot(1:lh, histf[1:lh]; legend=false)
        # title!("History objective value")
        xlabel!("Step")
        ylabel!("Objective Value "*L"f_{\beta}")
        annotate!(-0.1*length(histf[1:lh]),-(extrema(histf[1:lh])...)*0.95,
            (L"\times 10^3",10,:left))

        plt6 = plot(1:lh, histd[:, 1:lh]'/1e-3; legend=false)
        # title!("History distances")
        xlabel!("Step")
        ylabel!("Distances "*L"d_i"*" [mm]")

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







import Dragoon: SATrace

function analyse(hist,trace::Vector{SATrace},freqsplot;
        freqs=nothing,plotting=true,div=5,scale=1e9,ylim=[-0.05e4,3e4],
        freqticks=[21.95,22.0,22.05,22.1])

    tracex = hcat((x -> x.x).(trace)...)
    traced = hcat((x -> pos2dist(x.x)).(trace)...)

    tracexsol = hcat((x -> x.xsol).(trace)...)
    tracedsol = hcat((x -> pos2dist(x.xsol)).(trace)...)

    traceobj = (x -> x.obj).(trace)/1e3
    traceobjsol = (x -> x.objsol).(trace)/1e3

    tracet = (x -> x.τ).(trace)

    l = length(trace)
    n = length(tracex[:, 1])

    lh = length(hist[(x->x.objvalue).(hist).!=0.0])

    histx = hcat((x -> x.pos).(hist[lh:-1:1])...)
    histf = (x -> x.objvalue).(hist[lh:-1:1])/1e3
    histd = hcat((x -> pos2dist(x.pos)).(hist[lh:-1:1])...)

    if plotting
        plt1 = plot(freqsplot / scale, boost1d(pos2dist(tracexsol[:, 1]), freqsplot)/1e3;
            ylim=ylim/1e3, label="init", lc="blue", lw=2)

        if div != 0
            for i in 2:maximum([1, l ÷ div]):(l-1)
                plot!(freqsplot / scale, boost1d(pos2dist(tracexsol[:, i]), freqsplot)/1e3;
                    ylim=ylim/1e3, label="it.: " * string(i))
            end
        end

        plot!(freqsplot / scale, boost1d(pos2dist(tracexsol[:, l]), freqsplot)/1e3;
            ylim=ylim/1e3, label="final", lc="red", lw=2)

        if freqs !== nothing
            vline!([minimum(freqs), maximum(freqs)] / scale, c="black", linestyle=:dash,
                label="")
        end

        # title!("Boostfactor")
        xticks!(freqticks)
        xlabel!("Frequency [GHz]")
        ylabel!("Power Boost Factor "*L"\beta^2")
        annotate!([(minimum(freqsplot) / scale, 0.9 * ylim[2]/1e3,
            L"f_{\beta,F}:" * string(round(traceobjsol[l]*1e3, digits=0)), :left)])
        annotate!((minimum(freqsplot)+-(extrema(freqsplot)...)*0.1) / scale,
            -0.05*ylim[2]/1e3,(L"\times 10^3",10,:left))

        plt2 = plot(1:l, traceobjsol; legend=false)
        # title!("Solution objective trace")
        xlabel!("Iteration")
        ylabel!("Objective Value "*L"f_{\beta}")
        annotate!(-0.1*length(traceobjsol),-(extrema(traceobjsol)...)*0.95,
            (L"\times 10^3",10,:left))

        plt3 = plot(1:l, tracedsol'/1e-3; legend=false)
        # title!("Solution distance trace [m]")
        xlabel!("Iteration")
        ylabel!("Distances "*L"d_i"*" [mm]")

        plt4 = scatter(1:n, tracedsol[:, l]/1e-3; legend=false)
        # title!("Final distances")
        xlabel!("Disk")
        ylabel!("Distances "*L"d_i"*" [mm]")

        plt5 = plot(1:l, traceobj; legend=false)
        # title!("Thermal objective trace")
        xlabel!("Iteration")
        ylabel!("Objective value")

        plt6 = plot(1:lh, histd[:, 1:lh]'/1e-3; legend=false)
        # title!("Thermal distance trace")
        xlabel!("Iteration")
        ylabel!("Distances "*L"d_i"*" [mm]")

        plt7 = plot((x->x.objvalue).(hist[(x->x.objvalue).(hist).!=0.0])[end:-1:1][1:end]/1e3;
            legend=false)
        # title!("History")
        xlabel!("Step")
        ylabel!("Objective Value "*L"f_{\beta}")
        annotate!(-0.1*length(histf[1:lh]),-(extrema(histf[1:lh])...)*0.95,
            (L"\times 10^3",10,:left))

        display(plt1)
        display(plt2)
        display(plt3)
        display(plt4)
        display(plt5)
        display(plt6)
        display(plt7)

        return plt1, plt2, plt3, plt4, plt5, plt6, plt7
    else
        return tracex, traced, tracexsol, tracedsol, traceobj, traceobjsol, tracet
    end
end


import Dragoon.LSTrace

function analyse(hist, trace::Vector{LSTrace}, freqsplot;
        freqs=nothing, plotting=true, div=5, scale=1e9, ylim=[-0.05e4, 3e4],
        freqticks=[21.95,22.0,22.05,22.1])

    tracex = hcat((x -> x.x).(trace)...)
    traced = hcat((x -> pos2dist(x.x)).(trace)...)
    tracef = (x -> x.obj).(trace)/1e3
    traceg = hcat((x -> x.g).(trace)...)
    traceh = cat((x -> x.h).(trace)...; dims=3)

    l = length(trace)
    n = length(tracex[:, 1])

    lh = length(hist[(x->x.objvalue).(hist).!=0.0])

    histx = hcat((x -> x.pos).(hist[lh:-1:1])...)
    histf = (x -> x.objvalue).(hist[lh:-1:1])/1e3
    histd = hcat((x -> pos2dist(x.pos)).(hist[lh:-1:1])...)

    if plotting
        plt1 = plot(freqsplot / scale, boost1d(pos2dist(tracex[:, 1]), freqsplot)/1e3;
            ylim=ylim/1e3, label="init", lc="blue", lw=2)

        if div != 0
            for i in 2:maximum([1, l ÷ div]):(l-1)
                plot!(freqsplot / scale, boost1d(pos2dist(tracex[:, i]), freqsplot)/1e3;
                    ylim=ylim/1e3, label="it.: " * string(i))
            end
        end

        plot!(freqsplot / scale, boost1d(pos2dist(tracex[:, l]), freqsplot)/1e3;
            ylim=ylim/1e3, label="final", lc="red", lw=2)

        if freqs !== nothing
            vline!([minimum(freqs), maximum(freqs)] / scale, c="black", linestyle=:dash,
                label="")
        end
        # title!("Boostfactor")
        xticks!(freqticks)
        xlabel!("Frequency [GHz]")
        ylabel!("Power Boost Factor "*L"\beta^2")
        annotate!([(minimum(freqsplot) / scale, 0.9 * ylim[2]/1e3,
            L"f_{\beta,F}:" * string(round(tracef[l]*1e3, digits=0)), :left)])
        annotate!((minimum(freqsplot)+-(extrema(freqsplot)...)*0.1) / scale,
            -0.05*ylim[2]/1e3,(L"\times 10^3",10,:left))

        plt2 = plot(1:l, tracef; legend=false)
        # title!("Objective trace")
        xlabel!("Iteration")
        ylabel!("Objective Value "*L"f_{\beta}")
        annotate!(-0.1*length(tracef),-(extrema(tracef)...)*0.95,
            (L"\times 10^3",10,:left))

        plt3 = plot(1:l, traced'/1e-3; legend=false)
        # title!("Distance trace")
        xlabel!("Iteration")
        ylabel!("Distances "*L"d_i"*" [mm]")

        plt4 = scatter(1:n, traced[:, l]; legend=false)
        # title!("Final distances")
        xlabel!("Step")
        ylabel!("Distances "*L"d_i"*" [mm]")

        plt5 = plot(1:lh, histf[1:lh]; legend=false)
        # title!("History")
        xlabel!("Step")
        ylabel!("Objective Value "*L"f_{\beta}")
        annotate!(-0.1*length(histf[1:lh]),-(extrema(histf[1:lh])...)*0.95,
            (L"\times 10^3",10,:left))

        plt6 = plot(1:lh, histd[:, 1:lh]'/1e-3; legend=false)
        # title!("History distances")
        xlabel!("Step")
        ylabel!("Distances "*L"d_i"*" [mm]")

        display(plt1)
        display(plt2)
        display(plt3)
        display(plt4)
        display(plt5)
        display(plt6)

        return plt1, plt2, plt3, plt4, plt5, plt6
    else
        return tracex, traced, tracef, traceg, traceh
    end
end
