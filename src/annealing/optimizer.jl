
export simulatedAnnealing

"""
simulatedAnnealing(booster::Booster,hist::States,freqs::Array{Float64},
        rmax::Real,
        T::Callback,
        objFunction::Callback,
        unstuckinator::Callback;
        maxiter::Integer=Int(1e2),
        nreset::Int64=0,
        nresetterm::Int64=0,
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        traceevery::Int=1,
        resettimer::Bool=true,
        returntimes::Bool=false)
        
Perform simulated annealing optimization on a physical or analytical booster, returns
process trace.
[`Simulated Annealing`](https://en.wikipedia.org/wiki/Simulated_annealing)

# Arguments
- `booster::Booster`: Booster struct to operate on.
- `hist::Vector{State}`: External vector to save visited states.
- `freqs::Array{Float64}`: Frequencies to optimize for.
- `rmax::Real`: Maximum distance of neighbour to recent point.
- `T::Callback`: Temperature manager.
- `objFunction::Callback`: Objective function of the booster system.
- `unstuckinator::Callback`: Softresets optimization when stuck.
- `maxiter::Integer=Int(1e2)`: Maximum iterations (duh).
- `nreset::Int64=0`: Softreset to best found position after `nreset` bad neighbours.
    Set 0 to disable.
- `nresetterm::Int64=0`: Maximum allowed softresets before termination.
- `showtrace::Bool=false`: Display progress informations.
- `showevery::Integer=1`: Show trace every n iterations.
- `unstuckisiter::Bool=true`: Count iteration where unstucking occured towards maxiter.
- `traceevery::Int=1`: Write to trace every n iterations.
- `resettimer::Bool=true`: Reset timestamps of booster.
- `returntimes::Bool=false`: Return timestamps additionally to the trace.
"""
function simulatedAnnealing(booster::Booster,hist::States,freqs::Array{Float64},
        rmax::Real,
        T::Callback,
        objFunction::Callback,
        unstuckinator::Callback;
        maxiter::Integer=Int(1e2),
        nreset::Int64=0,
        nresetterm::Int64=0,
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        traceevery::Int=1,
        resettimer::Bool=true,
        returntimes::Bool=false)

    t0 = setTimes!(booster,resettimer)

    trace = Vector{SATrace}(undef,round(Int,maxiter/traceevery)+2)

    x = copy(booster.pos)
    objx = updateHist!(booster,hist,freqs,objFunction)
    
    xsol = copy(booster.pos)
    objsol = copy(objx)

    dx = Vector{Float64}(undef,length(x))

    iter = 1
    resetcounter = 0
    resetcounterterm = 0

    τ = copy(T.args[1])

    trace[1] = SATrace(x,objx,xsol,objsol,τ,iter,booster.timestamp,booster.summeddistance)

    while iter < maxiter
        τ = T.func(τ,T.args)

        findNeighbour!(dx,rmax)

        move(booster,x+dx; additive=false)
        obj_ = updateHist!(booster,hist,freqs,objFunction; force=true)

        if obj_ <= objx || rand() <= thermal(objx,obj_,τ)
            # x = copy(booster.pos)
            copyto!(x,booster.pos)
            objx = obj_
        end

        if obj_ <= objsol
            # xsol = copy(booster.pos)
            copyto!(xsol,booster.pos)
            objsol = obj_

            resetcounter = 0
            resetcounterterm = 0
        else
            resetcounter += 1
        end

        if nreset > 0 && resetcounter >= nreset
            showtrace && println("Resetting to best solution.")
            resetcounter = 0
            resetcounterterm += 1

            if nresetterm > 0 && resetcounterterm >= nresetterm
                showtrace && println("$nresetterm times resetted. Unstucking.")

                stuck = unstuckinator.func(booster, hist, freqs, objFunction,
                    objsol, unstuckinator.args; showtrace=showtrace)
    
                !unstuckisiter && (iter -= 1)
    
                stuck && break

                resetcounterterm = 0
            end
        end
    
        showtrace && iter%showevery == 0 && printSAIter(booster,objx,objsol,τ,iter)

        if Int(iter%traceevery)==0
            trace[Int(iter/traceevery)+1] = SATrace(x,objx,xsol,objsol,τ,iter,
                                    booster.timestamp,booster.summeddistance)
        end

        iter += 1
    end

    move(booster,xsol; additive=false)
    updateHist!(booster,hist,freqs,objFunction)
    
    idx = min(findlast(i->isassigned(trace,i),eachindex(trace))+1,length(trace))
    trace[idx] = SATrace(x,objx,xsol,objsol,τ,iter,booster.timestamp,booster.summeddistance)

    updateTimeStamp!(booster,:codetimestamp,resettimer,t0)

    term = printTermination(booster,hist,iter,maxiter,showtrace)

    return returntimes ? (trace[1:idx],term) : trace[1:idx]
end
