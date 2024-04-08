
export nelderMead

"""
    nelderMead(booster::Booster, hist::Vector{State}, freqs::Array{Float64},
        α::Float64, β::Float64, γ::Float64, δ::Float64, Δmin::Real,
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
        resettimer::Bool=true,
        returntimes::Bool=false)

Perform Nelder-Mead optimization on a physical or analytical booster, returns process trace.

# Arguments
- `booster::Booster`: Booster struct to operate on.
- `hist::Vector{State}`: External vector to save visited states.
- `freqs::Array{Float64}`: Frequencies to optimize for.
- `α::Float64`: Reflection parameter.
- `β::Float64`: Expansion parameter.
- `γ::Float64`: Contraction parameter.
- `δ::Float64`: Shrinking parameter.
- `Δmin::Real`: Minimum allowed simplex size.
- `objFunction::Callback`: Objective function of the booster system.
- `initSimplex`: Simplex initializer.
- `simplexObj`: Manager to iterate over all vertices when getting all objective values.
- `unstuckinator::Callback`: Softresets optimization when stuck.
- `maxiter::Integer=Int(1e2)`: Maximum iterations (duh).
- `showtrace::Bool=false`: Display progress informations.
- `showevery::Integer=1`: Show trace every n iterations.
- `unstuckisiter::Bool=true`: Count iteration where unstucking occured towards maxiter.
- `forcesimplexobj::Bool=false`: Forces optimizer to recalculate entire simplex in every
    iteration.
- `tracecentroid::Bool=false`: Additionally trace centroid of simplex.
- `resettimer::Bool=true`: Reset timestamps of booster.
- `returntimes::Bool=false`: Return timestamps additionally to the trace.
"""
function nelderMead(booster::Booster, hist::Vector{State}, freqs::Array{Float64},
        α::Float64, β::Float64, γ::Float64, δ::Float64, Δmin::Real,
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
        resettimer::Bool=true,
        returntimes::Bool=false)

    t0 = setTimes!(booster,resettimer)

    trace = Vector{NMTrace}(undef, floor(Int, maxiter / traceevery) + 1)

    x = initSimplex.func(booster.pos, initSimplex.args)
    f = simplexObj.func(x, collect(1:booster.ndisk+1), booster, hist, freqs,
        objFunction, simplexObj.args)

    i = 0
    while i < maxiter
        i += 1

        #sort
        sp = sortperm(f; rev=false)
        x[:, :] = x[:, sp]
        f[:] = f[sp]

        if Int(i % traceevery) == 0
            trace[Int(i / traceevery)] = NMTrace(copy(x), copy(f), zeros(booster.ndisk), 0.0,
                booster.timestamp, booster.summeddistance)
        end

        #centroid
        x_ = reshape(sum(x[:, 1:end-1]; dims=2), :) / booster.ndisk

        if tracecentroid
            move(booster, x_; additive=false)
            updateHist!(booster, hist, freqs, objFunction)

            trace[Int(i / traceevery)].x_ = x_
            trace[Int(i / traceevery)].obj_ = hist[1].objvalue
        end

        #reflection point, reflection
        xr = x_ + α * (x_ - x[:, end])
        move(booster, xr; additive=false)
        updateHist!(booster, hist, freqs, objFunction)

        fr = hist[1].objvalue

        if f[1] <= fr < f[end-1]
            x[:, end] = xr
            f[end] = fr
            #end

            #expansion
        elseif fr < f[1]
            xe = x_ + β * (xr - x_)

            move(booster, xe; additive=false)
            updateHist!(booster, hist, freqs, objFunction)

            fe = hist[1].objvalue

            if fe < fr
                x[:, end] = xe
                f[end] = fe
            else
                x[:, end] = xr
                f[end] = fr
            end
            #end

            #contraction
        else #if fr >= f[end-1]
            if f[end-1] <= fr < f[end] #outside contraction
                xoc = x_ + γ * (xr - x_)

                move(booster, xoc; additive=false)
                updateHist!(booster, hist, freqs, objFunction)

                foc = hist[1].objvalue

                if foc <= fr
                    x[:, end] = xoc
                    f[end] = foc
                else #shrink
                    for j in 2:booster.ndisk+1
                        v = x[:, 1] + δ * (x[:, j] - x[:, 1])
                        x[:, j] = v
                    end

                    f[2:end] = simplexObj.func(x, Vector(2:booster.ndisk+1),
                        booster, hist, freqs, objFunction, simplexObj.args)
                end
            else #inside contraction
                xic = x_ - γ * (xr - x_)

                move(booster, xic; additive=false)
                updateHist!(booster, hist, freqs, objFunction)

                fic = hist[1].objvalue

                if fic < fr
                    x[:, end] = xic
                    f[end] = fic
                else #shrink
                    for j in 2:booster.ndisk+1
                        v = x[:, 1] + δ * (x[:, j] - x[:, 1])
                        x[:, j] = v
                    end

                    f[2:end] = simplexObj.func(x, collect(2:booster.ndisk+1),
                        booster, hist, freqs, objFunction, simplexObj.args)
                end
            end
        end

        showtrace && i % showevery == 0 && printNMIter(booster, f, i)

        #iteration end

        if forcesimplexobj
            f[:] = simplexObj.func(x, collect(1:booster.ndisk+1), booster, hist, freqs,
                objFunction, simplexObj.args)
        end

        if getSimplexSize(x, f) < Δmin
            showtrace && println("Minimum simplex size reached.")

            stuck = unstuckinator.func(booster, hist, freqs, objFunction, simplexObj, x, f,
                unstuckinator.args; showtrace=showtrace)

            !unstuckisiter && (i -= 1)

            stuck && break
        end
    end

    #sort and trace the end result
    sp = sortperm(f; rev=false)
    x[:, :] = x[:, sp]
    f[:] = f[sp]

    trace[end] = NMTrace(x, f, zeros(booster.ndisk), 0.0, booster.timestamp,
        booster.summeddistance)

    x_ = reshape(sum(x[:, 1:end-1]; dims=2), :) / booster.ndisk

    if tracecentroid
        move(booster, x_; additive=false)
        updateHist!(booster, hist, freqs, objFunction)

        trace[end].x_ = x_
        trace[end].obj_ = hist[1].objvalue
    end

    #move to optimal point
    move(booster, x[:, argmin(f)]; additive=false)
    updateHist!(booster, hist, freqs, objFunction)

    updateTimeStamp!(booster,:codetimestamp,resettimer,t0)

    term = printTermination(booster, hist, i, maxiter, showtrace)

    return returntimes ? (trace[1:i+1], term) : trace[1:i+1]
end


function nelderMead_(booster::Booster, hist::Vector{State}, freqs::Array{Float64},
        α::Float64, β::Float64, γ::Float64, δ::Float64, Δmin::Real,
        objFunction::Callback,
        initSimplex::Callback,
        simplexObj::Callback,
        unstuckinator::Callback;
        maxiter::Integer=Int(1e2),
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        forcesimplexobj::Bool=false,
        traceevery::Int=1,
        resettimer::Bool=true,
        returntimes::Bool=false)

    t0 = setTimes!(booster, resettimer)

    trace = Vector{NMTrace}(undef, floor(Int, maxiter / traceevery) + 1)

    x = initSimplex.func(booster.pos, initSimplex.args)
    f = simplexObj.func(x, collect(1:booster.ndisk+1), booster, hist, freqs,
        objFunction, simplexObj.args)

    i = 0

    while i < maxiter
        i += 1

        sp = sortperm(f; rev=false) # sort
        x[:, :] = x[:, sp]
        f[:] = f[sp]

        x_ = reshape(sum(x[:, 1:end-1]; dims=2), :) / booster.ndisk
        xr = x_ + α * (x_ - x[:, end])
        xic = x_ - γ * (xr - x_)
        xoc = x_ + γ * (xr - x_)
        xe = x_ + β * (xr - x_)

        move(booster, xic; additive=false)
        fic = updateHist!(booster, hist, freqs, objFunction)

        move(booster, x_; additive=false)
        f_ = updateHist!(booster, hist, freqs, objFunction)

        move(booster,xoc; additive=false)
        foc = updateHist!(booster, hist, freqs, objFunction)

        move(booster,xoc; additive=false)
        fr = updateHist!(booster, hist, freqs, objFunction)

        if Int(i % traceevery) == 0
            trace[Int(i / traceevery)] = NMTrace(copy(x), copy(f), zeros(booster.ndisk), 0.0,
                booster.timestamp, booster.summeddistance)
            trace[Int(i / traceevery)].x_ = copy(x_)
            trace[Int(i / traceevery)].obj_ = f_
        end

        if f[1] <= fr < f[end-1]
            x[:, end] = xr
            f[end] = fr
        elseif fr < f[1]    # expansion
            move(booster, xe; additive=false)
            fe = updateHist!(booster, hist, freqs, objFunction)

            if fe < fr
                x[:, end] = xe
                f[end] = fe
            else
                x[:, end] = xr
                f[end] = fr
            end
        else # if fr >= f[end-1] # contraction
            if f[end-1] <= fr < f[end] # outside contraction
                if foc <= fr
                    x[:, end] = xoc
                    f[end] = foc
                else # shrink
                    for j in 2:booster.ndisk+1
                        v = x[:, 1] + δ * (x[:, j] - x[:, 1])
                        x[:, j] = v
                    end

                    f[2:end] = simplexObj.func(x, Vector(2:booster.ndisk+1),
                        booster, hist, freqs, objFunction, simplexObj.args)
                end
            else # inside contraction
                if fic < fr
                    x[:, end] = xic
                    f[end] = fic
                else # shrink
                    for j in 2:booster.ndisk+1
                        v = x[:, 1] + δ * (x[:, j] - x[:, 1])
                        x[:, j] = v
                    end

                    f[2:end] = simplexObj.func(x, collect(2:booster.ndisk+1),
                        booster, hist, freqs, objFunction, simplexObj.args)
                end
            end
        end

        showtrace && i % showevery == 0 && printNMIter(booster, f, i)   # iteration end

        if forcesimplexobj
            f[:] = simplexObj.func(x, collect(1:booster.ndisk+1), booster, hist, freqs,
                objFunction, simplexObj.args)
        end

        if getSimplexSize(x, f) < Δmin
            showtrace && println("Minimum simplex size reached.")

            stuck = unstuckinator.func(booster, hist, freqs, objFunction, simplexObj, x, f,
                unstuckinator.args; showtrace=showtrace)

            !unstuckisiter && (i -= 1)

            stuck && break
        end
    end

    sp = sortperm(f; rev=false) # sort and trace the end result
    x[:, :] = x[:, sp]
    f[:] = f[sp]

    trace[end] = NMTrace(x, f, zeros(booster.ndisk), 0.0, booster.timestamp,
        booster.summeddistance)

    trace[end].x_ = reshape(sum(x[:, 1:end-1]; dims=2), :) / booster.ndisk
    move(booster, trace[end].x_; additive=false)
    trace[end].obj_ = updateHist!(booster, hist, freqs, objFunction)
    
    move(booster, x[:, argmin(f)]; additive=false)  # move to optimal point
    updateHist!(booster, hist, freqs, objFunction)

    updateTimeStamp!(booster, :codetimestamp, resettimer, t0)

    term = printTermination(booster, hist, i, maxiter, showtrace)

    return returntimes ? (trace[1:i+1], term) : trace[1:i+1]
end


function nelderMeadLinesearch(booster::Booster, hist::Vector{State}, freqs::Array{Float64},
        α::Float64, β::Float64, γ::Float64, δ::Float64, Δmin::Real,αls::Real,
        objFunction::Callback,
        initSimplex::Callback,
        simplexObj::Callback,
        unstuckinator::Callback;
        maxiter::Integer=Int(1e2),
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        forcesimplexobj::Bool=false,
        traceevery::Int=1,
        resettimer::Bool=true,
        returntimes::Bool=false)

    t0 = setTimes!(booster, resettimer)

    trace = Vector{NMTrace}(undef, floor(Int, maxiter / traceevery) + 1)

    x = initSimplex.func(booster.pos, initSimplex.args)
    f = simplexObj.func(x, collect(1:booster.ndisk+1), booster, hist, freqs,
        objFunction, simplexObj.args)

    i = 0

    while i < maxiter
        i += 1

        sp = sortperm(f; rev=false) # sort
        x[:, :] = x[:, sp]
        f[:] = f[sp]

        if Int(i % traceevery) == 0
            trace[Int(i / traceevery)] = NMTrace(copy(x), copy(f), zeros(booster.ndisk), 0.0,
                booster.timestamp, booster.summeddistance)
        end

        x_ = reshape(sum(x[:, 1:end-1]; dims=2), :) / booster.ndisk
        xr_ = x_ + α * (x_ - x[:, end])
        xic_ = x_ - γ * (xr_ - x_)
        xoc_ = x_ + γ * (xr_ - x_)
        xe_ = x_ + β * (xr_ - x_)

        xic,fic = search(booster,hist,freqs,objFunction,x[:,end],xic_,αls)
        xr,fr = search(booster,hist,freqs,objFunction,xoc_,xr_,αls)

        if f[1] <= fr < f[end-1]
            x[:, end] = xr
            f[end] = fr
        elseif fr < f[1]    # expansion
            xe,fe = search(booster,hist,freqs,objFunction,xr_,xe_,αls)

            if fe < fr
                x[:, end] = xe
                f[end] = fe
            else
                x[:, end] = xr
                f[end] = fr
            end
        else # if fr >= f[end-1] # contraction
            if f[end-1] <= fr < f[end] # shrink
                for j in 2:booster.ndisk+1
                    v = x[:, 1] + δ * (x[:, j] - x[:, 1])
                    x[:, j] = v
                end

                f[2:end] = simplexObj.func(x, Vector(2:booster.ndisk+1),
                    booster, hist, freqs, objFunction, simplexObj.args)
            else # inside contraction
                if fic < fr
                    x[:, end] = xic
                    f[end] = fic
                else # shrink
                    for j in 2:booster.ndisk+1
                        v = x[:, 1] + δ * (x[:, j] - x[:, 1])
                        x[:, j] = v
                    end

                    f[2:end] = simplexObj.func(x, collect(2:booster.ndisk+1),
                        booster, hist, freqs, objFunction, simplexObj.args)
                end
            end
        end

        showtrace && i % showevery == 0 && printNMIter(booster, f, i)   # iteration end

        if forcesimplexobj
            f[:] = simplexObj.func(x, collect(1:booster.ndisk+1), booster, hist, freqs,
                objFunction, simplexObj.args)
        end

        if getSimplexSize(x, f) < Δmin
            showtrace && println("Minimum simplex size reached.")

            stuck = unstuckinator.func(booster, hist, freqs, objFunction, simplexObj, x, f,
                unstuckinator.args; showtrace=showtrace)

            !unstuckisiter && (i -= 1)

            stuck && break
        end
    end

    sp = sortperm(f; rev=false) # sort and trace the end result
    x[:, :] = x[:, sp]
    f[:] = f[sp]

    trace[end] = NMTrace(x, f, zeros(booster.ndisk), 0.0, booster.timestamp,
        booster.summeddistance)

    trace[end].x_ = reshape(sum(x[:, 1:end-1]; dims=2), :) / booster.ndisk
    move(booster, trace[end].x_; additive=false)
    trace[end].obj_ = updateHist!(booster, hist, freqs, objFunction)

    move(booster, x[:, argmin(f)]; additive=false)  # move to optimal point
    updateHist!(booster, hist, freqs, objFunction)

    updateTimeStamp!(booster, :codetimestamp, resettimer, t0)

    term = printTermination(booster, hist, i, maxiter, showtrace)

    return returntimes ? (trace[1:i+1], term) : trace[1:i+1]
end

function search(booster::Booster,hist::States,freqs::Vector{Float64},objFunction::Callback,
        x_start::Vector{Float64},x_stop::Vector{Float64},αls::Real)

    dx = x_stop-x_start; dx_l = pNorm(dx); dx /= dx_l
    nsteps = floor(dx_l/αls)

    move(booster,x_start; additive=false)
    updateHist!(booster,hist,freqs,objFunction)
    # f_start = updateHist!(booster,hist,freqs,objFunction)

    # f_best = Inf64; idx_best = 0

    # for i in nsteps:-1:1
    for i in 1:nsteps
        move(booster,αls*dx; additive=true)
        f_ = updateHist!(booster,hist,freqs,objFunction)

        # if f_ < f_best
        #     f_best = f_; idx_best = i
        # end
    end

    # x_best = copy(hist[idx_best].pos)
    
    move(booster,x_stop; additive=false)
    updateHist!(booster,hist,freqs,objFunction)
    # f_stop = updateHist!(booster,hist,freqs,objFunction)

    idx = argmin(x->x.objvalue,hist[1:nsteps+2])

    return hist[idx].pos, hist[idx].objvalue
end