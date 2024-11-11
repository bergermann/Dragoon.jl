
export adam



function adam(booster::Booster, hist::Vector{State}, freqs::Array{Float64},
        α::Float64,β1::Float64,β2::Float64,ϵ::Real,
        objFunction::Callback,
        derivator::Callback,
        unstuckinator::Callback;
        ϵgrad::Float64=0.0,
        maxiter::Integer=Int(1e2),
        showtrace::Bool=false,
        showevery::Integer=1,
        traceevery::Int=1,
        resettimer::Bool=true,
        returntimes::Bool=false)

    @assert (0 <= β1 < 1) && (0 <= β2 < 1) "β1, β2 need to be from the range [0,1)."

    t0 = setTimes!(booster,resettimer)

    trace = Vector{ATrace}(undef,round(Int,maxiter/traceevery)+2)

    m = zeros(booster.ndisk)
    v = zeros(booster.ndisk)

    g = zeros(booster.ndisk)
    h = nothing

    β1_ = β1; β2_ = β2; α_ = α

    θ = copy(booster.pos)

    iter = 0
    while iter < maxiter
        updateHist!(booster, hist, freqs, objFunction)

        derivator.func(g, h, booster, hist, freqs, objFunction, derivator.args)

        if iter%traceevery == 0
            trace[Int(iter/traceevery)+1] = ATrace(booster.pos,
                hist[1].objvalue,g,booster.timestamp,booster.summeddistance)
        end
        
        iter += 1

        showtrace && iter%showevery == 0 && println("Gradient norm: ",
                                                round(pNorm(g), sigdigits=3))

        #early stopping if descend is too slow
        pNorm(g) <= ϵgrad && ((showtrace && println("Gradient threshold reached.
                                                Terminating.")); break)

        m *= β1; m += (1-β1)*g
        v *= β2; @. v += (1-β2)*g^2

        β1_ *= β1; β2_ *= β2

        α_ = α*sqrt(1-β2_)/(1-β1_)
        @. θ -= α_*m/(sqrt(v)+ϵ)

        move(booster,θ; additive=false)

        showtrace && iter%showevery == 0 && printAIter(booster, hist, iter)

        # stuck = unstuckinator.func(booster, hist, freqs, objFunction,
        #     unstuckinator.args; showtrace=showtrace)
        
        # stuck && break
    end

    updateTimeStamp!(booster,:codetimestamp,resettimer,t0)

    idx = min(findlast(i->isassigned(trace,i),eachindex(trace))+1,length(trace))
    trace[idx] = ATrace(booster.pos,hist[1].objvalue,g,booster.timestamp,booster.summeddistance)

    term = printTermination(booster,hist,iter,maxiter,showtrace)

    return returntimes ? (trace[1:idx], term) : trace[1:idx]
end


