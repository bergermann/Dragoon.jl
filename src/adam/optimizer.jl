



function adam(booster::Booster, hist::Vector{State}, freqs::Array{Float64},
        α::Float64,β1::Float64,β2::Float64,ϵ::Real,
        objFunction::Callback,
        derivator::Callback,
        unstuckinator::Callback;
        ϵgrad::Float64=0.0,
        maxiter::Integer=Int(1e2),
        showtrace::Bool=false,
        showevery::Integer=1,
        unstuckisiter::Bool=true,
        resettimer::Bool=true,
        returntimes::Bool=false)

    @assert (0 <= β1 < 1) && (0 <= β2 < 1) "β1, β2 need to be from the range [0,1)."

    t0 = setTimes!(booster,resettimer)

    trace = Vector{ATrace}(undef,maxiter+1)

    m = zeros(booster.ndisk)
    v = zeros(booster.ndisk)

    g = zeros(booster.ndisk)
    h = nothing

    β1_ = β1; β2_ = β2; α_ = α

    θ = copy(booster.pos)

    i = 0
    while i < maxiter
        i += 1

        updateHist!(booster, hist, freqs, objFunction)

        derivator.func(g, h, booster, hist, freqs, objFunction, derivator.args)

        trace[i] = ATrace(booster.pos,hist[1].objvalue,copy(g),
            booster.timestamp,booster.summeddistance)

        showtrace && i%showevery == 0 && println("Gradient norm: ",
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

        # showtrace && i%showevery == 0 && printIter(booster, hist, i, k)

        stuck = unstuckinator.func(booster, hist, freqs, objFunction,
            unstuckinator.args; showtrace=showtrace)

        !unstuckisiter && (i -= 1) #reset iteration count if false

        stuck && break

        #end of iteration
    end

    updateTimeStamp!(booster,:codetimestamp,resettimer,t0)

    idx = min(findlast(i->isassigned(trace,i),eachindex(trace))+1,length(trace))
    trace[idx] = LATrace(booster.pos,hist[1].objvalue,g,h,booster.timestamp, booster.summeddistance)

    term = printTermination(booster,hist,i,maxiter,showtrace)

    return returntimes ? (trace[1:idx], term) : trace[1:idx]
end


