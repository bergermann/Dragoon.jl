###     main file, optimization algorithms

export linesearch, nelderMead

import Dates: now, UTC

function linesearch(booster::AnalyticalBooster,hist::Vector{State},freqs::Array{Float64},
                    α::Float64,
                    objFunction::Tuple{Function,Vector},
                    solver::Tuple{Function,Vector},
                    derivatives::Tuple{Function,Vector},
                    step::Tuple{Function,Vector},
                    search::Tuple{Function,Vector},
                    unstuckinator::Tuple{Function,Vector};
                    ϵgrad::Float64=0,maxiter::Integer=Int(1e2),
                    showtrace::Bool=false,showevery::Integer=1,
                    unstuckisiter::Bool=true)

    t = now(UTC)

    trace = Vector{LSTrace}(undef,maxiter+1)

    p = zeros(booster.ndisk)
    g = zeros(booster.ndisk)
    h = zeros(booster.ndisk,booster.ndisk)

    i = 0
    while i < maxiter
        i += 1

        #calculate derivatives and step direction
        updateHist!(booster,hist,freqs,objFunction)

        derivatives[1](g,h,booster,hist,freqs,objFunction,derivatives[2]...)

        trace[i] = LSTrace(booster.pos,hist[1].objvalue,copy(g),copy(h),
                                    booster.timestamp,booster.summedtraveltime)

        solver[1](p,g,h,trace,i,solver[2]...)

        showtrace && println("Gradient norm: ",round(pNorm(g),sigdigits=3))

        #early stopping if descend is too slow
        pNorm(g) <= ϵgrad && (println("Gradient threshold reached.
                                                        Terminating."); break)

        #determine steplength and/or normalize p
        updateHist!(booster,hist,freqs,objFunction)

        step[1](p,α,booster,hist,freqs,objFunction,step[2]...;
                                                        showtrace=showtrace)

        #perform linesearch
        updateHist!(booster,hist,freqs,objFunction)

        k = search[1](p,α,booster,hist,freqs,objFunction,search[2]...;
                                                        showtrace=showtrace)

        showtrace && i%showevery == 0 && printIter(booster,hist,i,k)

        #perform unstucking
        if k == 0
            stuck = unstuckinator[1](booster,hist,freqs,objFunction,
                                        unstuckinator[2]...; showtrace=showtrace)

            !unstuckisiter && (i -= 1) #reset iteration count if false

            stuck && break
        end

        #end of iteration
    end

    booster.codetimestamp = canonicalize(now(UTC)-t)
    printTermination(booster,hist,i,maxiter)

    trace[i+1] = LSTrace(booster.pos,hist[1].objvalue,g,h,
                                booster.timestamp,booster.summedtraveltime)

    return trace[1:i+1]
end

function nelderMead(booster::AnalyticalBooster,hist::Vector{State},freqs::Array{Float64},
                    α::Float64,β::Float64,γ::Float64,δ::Float64,
                    objFunction::Tuple{Function,Vector},
                    initSimplex::Tuple{Function,Vector},
                    simplexObj::Tuple{Function,Vector},
                    unstuckinator::Tuple{Function,Vector};
                    maxiter::Integer=Int(1e2),
                    showtrace::Bool=false,
                    showevery::Integer=1,
                    unstuckisiter::Bool=true,
                    tracecentroid::Bool=false,
                    traceevery::Int=1)

    t = now(UTC)

    trace = Vector{NMTrace}(undef,Int(maxiter/traceevery)+1)

    x = zeros(booster.ndisk,booster.ndisk+1)
    f = zeros(booster.ndisk+1)

    initSimplex[1](booster.pos,x,initSimplex[2]...)
    f = simplexObj[1](x,Vector(1:booster.ndisk+1),booster,hist,freqs,
                                            objFunction,simplexObj[2]...)

    i = 0
    while i < maxiter
        i += 1

        #sort
        sp = sortperm(f; rev=false)
        x = x[:,sp]
        f = f[sp]

        if Int(i%traceevery)==0
            trace[Int(i/traceevery)] = NMTrace(x,f,zeros(booster.ndisk),0.,
                                    booster.timestamp,booster.summedtraveltime)
        end

        #centroid
        x_ = reshape(sum(x[:,1:end-1]; dims=2),:)/booster.ndisk
        if tracecentroid
            move(booster,x_; additive=false)
            updateHist!(booster,hist,freqs,objFunction)

            trace[i].x_ = x_
            trace[i].obj = hist[1].objvalue
        end

        #reflection point, reflection
        xr = x_+α*(x_-x[:,end])
        move(booster,xr; additive=false)
        updateHist!(booster,hist,freqs,objFunction)
        fr = hist[1].objvalue

        if f[1] <= fr < f[end-1]
            x[:,end] = xr
            f[end] = fr
        #end

        #expansion
        elseif fr < f[1]
            xe = x_+β*(xr-x_)

            move(booster,xe; additive=false)
            updateHist!(booster,hist,freqs,objFunction)

            fe = hist[1].objvalue

            if fe < fr
                x[:,end] = xe
                f[end] = fe
            else
                x[:,end] = xr
                f[end] = fr
            end
        #end

        #contraction
        else #if fr >= f[end-1]
            if f[end-1] <= fr < f[end] #outside contraction
                xoc = x_+γ*(xr-x_)

                move(booster,xoc; additive=false)
                updateHist!(booster,hist,freqs,objFunction)

                foc = hist[1].objvalue

                if foc <= fr
                    x[:,end] = xoc
                    f[end] = foc
                else #shrink
                    for j in 2:booster.ndisk+1
                        v = x[:,1]+δ*(x[:,j]-x[:,1])
                        x[:,j] = v
                    end
                    f[2:end] = simplexObj[1](x,Vector(2:booster.ndisk+1),
                                booster,hist,freqs,objFunction,simplexObj[2]...)
                end
            else #inside contraction
                xic = x_-γ*(xr-x_)

                move(booster,xic; additive=false)
                updateHist!(booster,hist,freqs,objFunction)

                fic = hist[1].objvalue

                if fic < fr
                    x[:,end] = xic
                    f[end] = fic
                else #shrink
                    for j in 2:booster.ndisk+1
                        v = x[:,1]+δ*(x[:,j]-x[:,1])
                        x[:,j] = v
                    end
                    f[2:end] = simplexObj[1](x,Vector(2:booster.ndisk+1),
                                booster,hist,freqs,objFunction,simplexObj[2]...)
                end
            end
        end

        showtrace && i%showevery==0 && printNMIter(booster,f,i)

        #iteration end
    end

    #sort and trace the end result
    sp = sortperm(f; rev=false)
    x = x[:,sp]
    f = f[sp]

    trace[Int(i/traceevery)+1] = NMTrace(x,f,zeros(booster.ndisk),0.,booster.timestamp,
                                                booster.summedtraveltime)

    x_ = reshape(sum(x[:,1:end-1]; dims=2),:)/booster.ndisk
    if tracecentroid
        move(booster,x_; additive=false)
        updateHist!(booster,hist,freqs,objFunction)

        trace[Int(i/traceevery)+1].x_ = x_
        trace[Int(i/traceevery)+1].obj = hist[1].objvalue
    end

    #move to optimal point
    move(booster,x[:,argmin(f)]; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    booster.codetimestamp = canonicalize(now(UTC)-t)
    printTermination(booster,hist,i,maxiter)

    return trace[1:Int(i/traceevery)+1]
end


function nelderMead(booster::PhysicalBooster,hist::Vector{State},freqs::Array{Float64},
        α::Float64,β::Float64,γ::Float64,δ::Float64,
        objFunction::Tuple{Function,Vector},
        initSimplex::Tuple{Function,Vector},
        simplexObj::Tuple{Function,Vector},
        unstuckinator::Tuple{Function,Vector};
        maxiter::Integer=Int(1e2),
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        tracecentroid::Bool=false,
        traceevery::Int=1)

    t = now(UTC)

    trace = Vector{NMTrace}(undef,Int(maxiter/traceevery)+1)

    x = zeros(booster.ndisk,booster.ndisk+1)
    f = zeros(booster.ndisk+1)

    initSimplex[1](booster.pos,x,initSimplex[2]...)
    f = simplexObj[1](x,Vector(1:booster.ndisk+1),booster,hist,freqs,
                                objFunction,simplexObj[2]...)

    i = 0
    while i < maxiter
        i += 1

        #sort
        sp = sortperm(f; rev=false)
        x = x[:,sp]
        f = f[sp]

        if Int(i%traceevery)==0
            trace[Int(i/traceevery)] = NMTrace(x,f,zeros(booster.ndisk),0.,
                            booster.timestamp,booster.summedtraveltime)
        end

        #centroid
        x_ = sum(x[:,1:end-1]; dims=2)/booster.ndisk
        if tracecentroid
            move(booster,x_; additive=false)
            updateHist!(booster,hist,freqs,objFunction)

            trace[i].x_ = x_
            trace[i].obj = hist[1].objvalue
        end

        #reflection point, reflection
        xr = x_+α*(x_-x[:,end])
        move(booster,xr; additive=false)
        updateHist!(booster,hist,freqs,objFunction)
        fr = hist[1].objvalue

        if f[1] <= fr < f[end-1]
            x[:,end] = xr
            f[end] = fr
        elseif fr < f[1]    #expansion
            xe = x_+β*(xr-x_)

            move(booster,xe; additive=false)
            updateHist!(booster,hist,freqs,objFunction)

            fe = hist[1].objvalue

            if fe < fr
                x[:,end] = xe
                f[end] = fe
            else
                x[:,end] = xr
                f[end] = fr
            end
        else #if fr >= f[end-1] #contraction
            if f[end-1] <= fr < f[end] #outside contraction
                xoc = x_+γ*(xr-x_)

                move(booster,xoc; additive=false)
                updateHist!(booster,hist,freqs,objFunction)

                foc = hist[1].objvalue

                if foc <= fr
                    x[:,end] = xoc
                    f[end] = foc
                else #shrink
                    for j in 2:booster.ndisk+1
                        v = x[:,1]+δ*(x[:,j]-x[:,1])
                        x[:,j] = v
                    end
                    f[2:end] = simplexObj[1](x,Vector(2:booster.ndisk+1),
                                booster,hist,freqs,objFunction,simplexObj[2]...)
                end
            else #inside contraction
                xic = x_-γ*(xr-x_)

                move(booster,xic; additive=false)
                updateHist!(booster,hist,freqs,objFunction)

                fic = hist[1].objvalue

                if fic < fr
                    x[:,end] = xic
                    f[end] = fic
                else #shrink
                    for j in 2:booster.ndisk+1
                        v = x[:,1]+δ*(x[:,j]-x[:,1])
                        x[:,j] = v
                    end
                    f[2:end] = simplexObj[1](x,Vector(2:booster.ndisk+1),
                                booster,hist,freqs,objFunction,simplexObj[2]...)
                end
            end
        end

        showtrace && i%showevery==0 && printNMIter(booster,f,i)

        #iteration end
    end

    #sort and trace the end result
    sp = sortperm(f; rev=false)
    x = x[:,sp]
    f = f[sp]

    trace[Int(i/traceevery)+1] = NMTrace(x,f,zeros(booster.ndisk),0.,booster.timestamp,
                                    booster.summedtraveltime)

    x_ = sum(x[:,1:end-1]; dims=2)/booster.ndisk
    if tracecentroid
        move(booster,x_; additive=false)
        updateHist!(booster,hist,freqs,objFunction)

        trace[Int(i/traceevery)+1].x_ = x_
        trace[Int(i/traceevery)+1].obj = hist[1].objvalue
    end

    #move to optimal point
    move(booster,x[:,argmin(f)]; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    booster.codetimestamp = canonicalize(now(UTC)-t)
    printTermination(booster,hist,i,maxiter)

    return trace[1:Int(i/traceevery)+1]
end
