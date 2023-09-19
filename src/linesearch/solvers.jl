### solvers to determine directions

export solverSteepestDescent, solverNewton, solverBFGS,
    solverHybrid
export SolverSteep, SolverNewton, SolverBFGS, SolverHybrid

import LinearAlgebra: cholesky


# args = ()
function solverSteepestDescent(booster::Booster,hist::Vector{State},
        freqs::Vector{Float64},objFunction::Callback,p,g,h,trace,i,args)
    
    p[:] = -g
end

const SolverSteep = Callback(solverSteepestDescent)



# args = (mode,)
function solverNewton(booster::Booster,hist::Vector{State},freqs::Vector{Float64},
        objFunction::Callback,
        p::Vector{Float64},g::Vector{Float64},h::Matrix{Float64},
        trace::Vector{LSTrace},i::Int,args::Tuple{String})
    
    if args[1] == "cholesky"
        try
            C = cholesky(h)

            p[:] = -inv(C.U)*inv(C.L)*g
        catch e
            println("Hessian not Cholesky decomposable: ", e)
            println("Falling back to standard inversion.")
            
            p[:] = -inv(h)*g
        end
    else
        p[:] = -inv(h)*g
    end
end

const SolverNewton(mode::String) = Callback(solverNewton,(mode,))



# args = (h0,)
function solverBFGS(booster::Booster,hist::Vector{State},freqs::Vector{Float64},
        objFunction::Callback,p,g,h,trace,i,args)
    if i == 1
        trace[i].h = args[1]

        p[:] = inv(h0)*g

        return
    end

    y = g - trace[i-1].g
    s = trace[i].x - trace[i-1].x
    h_ = trace[i-1].h
    h = h_ + y*y'/(y'*s) - h_*s*s'*h_'/(s'*h_*s)
    trace[i].h = h

    p[:] = inv(h)*g
end

const SolverBFGS(h0) = Callback(solverBFGS,(h0,))




# args = (mode,ϵls,αtest,ntest)
function solverHybrid(booster::Booster,hist::Vector{State},freqs::Vector{Float64},
        objFunction::Callback,p::Vector{Float64},g::Vector{Float64},h::Matrix{Float64},
        trace::Vector{Dragoon.LSTrace},i::Int,
        (mode,ϵls,αtest,ntest)::Tuple{String,Real,Real,Int})

    if mode == "cholesky"
        try
            C = cholesky(h)

            p[:] = inv(C.U)*inv(C.L)*g
        catch e
            println("Hessian not Cholesky decomposable: ", e)
            println("Falling back to standard inversion.")
            
            p[:] = inv(h)*g
        end
    else
        p[:] = inv(h)*g
    end

    p[:] = p/pNorm(p)

    p0 = copy(booster.pos)
    updateHist!(booster,hist,freqs,objFunction)

    # test forward
    for i in ntest
        move(booster,αtest*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue > hist[i+1].objvalue+ϵls
            move(booster,p0)
            updateHist!(booster,hist,freqs,objFunction)

            break
        end

        if i == ntest
            move(booster,p0)
            updateHist!(booster,hist,freqs,objFunction)

            return
        end
    end

    # test backward
    for i in ntest
        move(booster,-αtest*p; additive=true)
        updateHist!(booster,hist,freqs,objFunction)

        if hist[1].objvalue > hist[i+1].objvalue+ϵls
            move(booster,p0)
            updateHist!(booster,hist,freqs,objFunction)

            break
        end

        if i == ntest
            move(booster,p0)
            updateHist!(booster,hist,freqs,objFunction)

            p[:] = -p[:]
            
            return
        end
    end

    # fall back to steepest descend
    p[:] = -g

    return
end

const SolverHybrid(mode::String,ϵls::Real,αtest::Real,ntest::Int) =
    Callback(solverHybrid,(mode,ϵls,αtest,ntest))
