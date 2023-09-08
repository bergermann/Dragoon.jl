###     main file, optimization algorithms

export linesearch, nelderMead, simulatedAnnealing

function linesearch(booster::Booster,hist::Vector{State},freqs::Array{Float64},
                    α::Float64,
                    objFunction::Callback,
                    solver::Callback,
                    derivator::Callback,
                    step::Callback,
                    search::Callback,
                    unstuckinator::Callback;
                    ϵgrad::Float64=0.,
                    maxiter::Integer=Int(1e2),
                    showtrace::Bool=false,
                    showevery::Integer=1,
                    unstuckisiter::Bool=true,
                    resettimer::Bool=true)

    if hasproperty(booster,:startingtime) && resettimer
        println("Resetting starting time.")
        
        booster.startingtime = unow()
        booster.timestamp = unow()
    elseif hasproperty(booster,:codetimestamp)
        t0 = unow()

        if resettimer
            booster.timestamp = DateTime(0)
        end
    end

    trace = Vector{LSTrace}(undef,maxiter+1)

    p = zeros(booster.ndisk)
    g = zeros(booster.ndisk)
    h = zeros(booster.ndisk,booster.ndisk)

    i = 0
    while i < maxiter
        i += 1

        #calculate derivatives and step direction
        updateHist!(booster,hist,freqs,objFunction)

        derivator.func(g,h,booster,hist,freqs,objFunction,derivator.args)

        # display(g); display(h)

        trace[i] = LSTrace(booster.pos,hist[1].objvalue,copy(g),copy(h),
                                    booster.timestamp,booster.summeddistance)

        solver.func(p,g,h,trace,i,solver.args)

        showtrace && i%showevery == 0 && println("Gradient norm: ",round(pNorm(g),sigdigits=3))

        #early stopping if descend is too slow
        pNorm(g) <= ϵgrad && (println("Gradient threshold reached.
                                                        Terminating."); break)

        #determine steplength and/or normalize p
        updateHist!(booster,hist,freqs,objFunction)

        step.func(p,α,booster,hist,freqs,objFunction,step.args; showtrace=showtrace)

        #perform linesearch
        updateHist!(booster,hist,freqs,objFunction)

        k = search.func(p,α,booster,hist,freqs,objFunction,search.args; showtrace=showtrace)

        showtrace && i%showevery == 0 && printIter(booster,hist,i,k)

        #perform unstucking
        if k == 0
            stuck = unstuckinator.func(booster,hist,freqs,objFunction,
                                        unstuckinator.args; showtrace=showtrace)

            !unstuckisiter && (i -= 1) #reset iteration count if false

            stuck && break
        end

        #end of iteration
    end


    if hasproperty(booster,:codetimestamp)
        if resettimer
            booster.codetimestamp = DateTime(0)
        end

        booster.codetimestamp += unow()-t0
    end

    printTermination(booster,hist,i,maxiter)

    trace[i+1] = LSTrace(booster.pos,hist[1].objvalue,g,h,
                                booster.timestamp,booster.summeddistance)

    return trace[1:i+1]
end





function nelderMead(booster::Booster,hist::Vector{State},freqs::Array{Float64},
                    α::Float64,β::Float64,γ::Float64,δ::Float64,
                    objFunction::Callback,
                    initSimplex::Callback,
                    simplexObj::Callback,
                    unstuckinator::Callback;
                    maxiter::Integer=Int(1e2),
                    showtrace::Bool=false,
                    showevery::Integer=1,
                    unstuckisiter::Bool=true,
                    forcesimplexobj::Bool=false,
                    tracecentroid::Bool=false,
                    traceevery::Int=1,
                    resettimer::Bool=true)

    if hasproperty(booster,:startingtime) && resettimer
        println("Resetting starting time.")
        
        booster.startingtime = unow()
        booster.timestamp = unow()
    elseif hasproperty(booster,:codetimestamp)
        t0 = unow()

        if resettimer
            booster.timestamp = DateTime(0)
        end
    end

    trace = Vector{NMTrace}(undef,floor(Int,maxiter/traceevery)+1)

    # x = zeros(booster.ndisk,booster.ndisk+1)
    # f = zeros(booster.ndisk+1)

    x = initSimplex.func(booster.pos,initSimplex.args)
    f = simplexObj.func(x,collect(1:booster.ndisk+1),booster,hist,freqs,
                                            objFunction,simplexObj.args)

    i = 0
    while i < maxiter
        i += 1

        #sort
        sp = sortperm(f; rev=false)
        x[:,:] = x[:,sp]
        f[:] = f[sp]

        if Int(i%traceevery)==0
            trace[Int(i/traceevery)] = NMTrace(x,f,zeros(booster.ndisk),0.,
                                    booster.timestamp,booster.summeddistance)
        end

        #centroid
        x_ = reshape(sum(x[:,1:end-1]; dims=2),:)/booster.ndisk
        if tracecentroid
            move(booster,x_; additive=false)
            updateHist!(booster,hist,freqs,objFunction)

            trace[i].x_ = x_
            trace[i].obj_ = hist[1].objvalue
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

                    f[2:end] = simplexObj.func(x,Vector(2:booster.ndisk+1),
                                booster,hist,freqs,objFunction,simplexObj.args)
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

                    f[2:end] = simplexObj.func(x,collect(2:booster.ndisk+1),
                                booster,hist,freqs,objFunction,simplexObj.args)
                end
            end
        end

        showtrace && i%showevery==0 && printNMIter(booster,f,i)

        #iteration end

        if forcesimplexobj
            f[:] = simplexObj.func(x,collect(1:booster.ndisk+1),booster,hist,freqs,
                objFunction,simplexObj.args)
        end
    end

    #sort and trace the end result
    sp = sortperm(f; rev=false)
    x[:,:] = x[:,sp]
    f[:] = f[sp]

    trace[end] = NMTrace(x,f,zeros(booster.ndisk),0.,booster.timestamp,
                                                booster.summeddistance)

    x_ = reshape(sum(x[:,1:end-1]; dims=2),:)/booster.ndisk

    if tracecentroid
        move(booster,x_; additive=false)
        updateHist!(booster,hist,freqs,objFunction)

        trace[end].x_ = x_
        trace[end].obj_ = hist[1].objvalue
    end

    #move to optimal point
    move(booster,x[:,argmin(f)]; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    if hasproperty(booster,:codetimestamp)
        if resettimer
            booster.codetimestamp = DateTime(0)
        end

        booster.codetimestamp += unow()-t0
    end

    printTermination(booster,hist,i,maxiter)

    return trace[1:Int(i/traceevery)+1]
end





function simulatedAnnealing(booster::Booster,hist::Vector{State},freqs::Array{Float64},
        τ::Array{Float64},rmax::Float64,
        objFunction::Callback,
        unstuckinator::Callback;
        maxiter::Integer=Int(1e2),
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        traceevery::Int=1,
        resettimer::Bool=true)

    if hasproperty(booster,:startingtime) && resettimer
        println("Resetting starting time.")
        
        booster.startingtime = unow()
        booster.timestamp = unow()
    elseif hasproperty(booster,:codetimestamp)
        t0 = unow()

        if resettimer
            booster.timestamp = DateTime(0)
        end
    end

    trace = Vector{SATrace}(undef,round(Int,maxiter/traceevery)+1)

    updateHist!(booster,hist,freqs,objFunction)
    
    x = copy(booster.pos)
    objx = hist[1].objvalue
    
    xsol = copy(booster.pos)
    objsol = hist[1].objvalue

    iter = 0
    i = 0
    n_τ = length(τ)

    while iter < maxiter && i < n_τ
        if Int(iter%traceevery)==0
            trace[Int(iter/traceevery)+1] = SATrace(x,objx,xsol,objsol,τ[i+1],iter,
                                    booster.timestamp,booster.summeddistance)
        end

        iter += 1
        i += 1

        # updateHist!(booster,hist,freqs,objFunction; force=true)
        move(booster,x+findNeighbour(booster,rmax); additive=false)
        updateHist!(booster,hist,freqs,objFunction; force=true)

        if hist[1].objvalue <= objx || rand() <= thermal(objx,hist[1].objvalue,τ[i])
            x = copy(booster.pos)
            objx = hist[1].objvalue
        end

        if hist[1].objvalue <= objsol
            xsol = copy(booster.pos)
            objsol = hist[1].objvalue
        end

        showtrace && i%showevery == 0 && printSAIter(booster,objx,objsol,τ[i],iter)

        ## unstucking, alter i
    end



    move(booster,xsol; additive=false)
    updateHist!(booster,hist,freqs,objFunction)

    if hasproperty(booster,:codetimestamp)
        if resettimer
            booster.codetimestamp = DateTime(0)
        end

        booster.codetimestamp += unow()-t0
    end

    printTermination(booster,hist,i,maxiter)

    return trace[1:min(round(Int,iter/traceevery)+1,length(trace))]
    # return trace
end
