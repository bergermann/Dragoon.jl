
export linesearch

"""
    linesearch(booster::Booster, hist::Vector{State}, freqs::Array{Float64},
        α::Float64,
        objFunction::Callback,
        solver::Callback,
        derivator::Callback,
        step::Callback,
        search::Callback,
        unstuckinator::Callback;
        ϵgrad::Float64=0.0,
        maxiter::Integer=Int(1e2),
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        resettimer::Bool=true,
        returntimes::Bool=false)

Perform linesearch optimization on a physical or analytical booster, returns process trace.
[`Linesearch`](https://en.wikipedia.org/wiki/Line_search)

# Arguments
- `booster::Booster`: Booster struct to operate on.
- `hist::Vector{State}`: External vector to save visited states.
- `freqs::Array{Float64}`: Frequencies to optimize for.
- `α::Float64`: Step length factor of linesearch.
- `objFunction::Callback`: Objective function of the booster system.
- `solver::Callback`: Solver that decides the search direction, see linesearch/solvers.jl.
- `derivator::Callback`: Calculates first/second derivatives of objective function.
- `step::Callback`: Optional normalization of step vector.
- `search::Callback`: Linesearch strategy.
- `unstuckinator::Callback`: Softresets optimization when stuck.
- `ϵgrad::Float64=0.0`: Minimum required slope of descend direction.
- `maxiter::Integer=Int(1e2)`: Maximum iterations (duh).
- `showtrace::Bool=false`: Display progress informations.
- `showevery::Integer=1`: Show trace every n iterations.
- `unstuckisiter::Bool=true`: Count iteration where unstucking occured towards maxiter.
- `resettimer::Bool=true`: Reset timestamps of booster.
- `returntimes::Bool=false`: Return timestamps additionally to the trace.
"""
function linesearch(booster::Booster, hist::Vector{State}, freqs::Array{Float64},
        α::Float64,
        objFunction::Callback,
        solver::Callback,
        derivator::Callback,
        step::Callback,
        search::Callback,
        unstuckinator::Callback;
        ϵgrad::Float64=0.0,
        maxiter::Integer=Int(1e2),
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        resettimer::Bool=true,
        returntimes::Bool=false)

    t0 = setTimes!(booster,resettimer)

    trace = Vector{LSTrace}(undef,maxiter+1)

    p = zeros(booster.ndisk)
    g = zeros(booster.ndisk)
    h = zeros(booster.ndisk,booster.ndisk)

    i = 0
    while i < maxiter
        i += 1

        #calculate derivatives and step direction
        updateHist!(booster, hist, freqs, objFunction)

        derivator.func(g, h, booster, hist, freqs, objFunction, derivator.args)

        trace[i] = LSTrace(booster.pos, hist[1].objvalue, copy(g), copy(h),
            booster.timestamp, booster.summeddistance)

        solver.func(booster, hist, freqs, objFunction, p, g, h, trace, i, solver.args)

        showtrace && i%showevery == 0 && println("Gradient norm: ",
                                                round(norm(g), sigdigits=3))

        #early stopping if descend is too slow
        norm(g) <= ϵgrad && ((showtrace && println("Gradient threshold reached.
                                                Terminating.")); break)

        #determine steplength and/or normalize p
        updateHist!(booster, hist, freqs, objFunction)

        step.func(p, α, booster, hist, freqs, objFunction, step.args; showtrace=showtrace)

        #perform linesearch
        updateHist!(booster, hist, freqs, objFunction)

        k = search.func(p, α, booster, hist, freqs, objFunction, search.args;
            showtrace=showtrace)

        showtrace && i%showevery == 0 && printIter(booster, hist, i, k)

        #perform unstucking
        if k == 0
            stuck = unstuckinator.func(booster, hist, freqs, objFunction,
                unstuckinator.args; showtrace=showtrace)

            !unstuckisiter && (i -= 1) #reset iteration count if false

            stuck && break
        end

        #end of iteration
    end

    updateTimeStamp!(booster,:codetimestamp,resettimer,t0)

    idx = min(findlast(i->isassigned(trace,i),eachindex(trace))+1,length(trace))
    trace[idx] = LSTrace(booster.pos,hist[1].objvalue,g,h,booster.timestamp, booster.summeddistance)

    term = printTermination(booster,hist,i,maxiter,showtrace)

    return returntimes ? (trace[1:idx], term) : trace[1:idx]
end


