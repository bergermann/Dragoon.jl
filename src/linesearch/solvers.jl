### solvers to determine directions

export solverSteepestDescent, solverNewton, solverBFGS,
    solverHybrid
export SolverSteep, SolverNewton, SolverBFGS, SolverHybrid

import LinearAlgebra: cholesky



"""
    solverSteepestDescent(booster::Booster,hist::Vector{State},
        freqs::Vector{Float64},objFunction::Callback,p,g,h,trace,i,args)

Return direction of steepest descent, ``-g``.
"""
function solverSteepestDescent(booster::Booster,hist::Vector{State},
        freqs::Vector{Float64},objFunction::Callback,p,g,h,trace,i,args)
    
    p[:] = -g
end

"""
    SolverSteep(mode)

Callback for steepest descend solver in [`linesearch`](@ref). See [`solverSteep`](@ref).
"""
SolverSteep = Callback(solverSteepestDescent)


"""
    function solverNewton(booster::Booster,hist::Vector{State},freqs::Vector{Float64},
        objFunction::Callback,
        p::Vector{Float64},g::Vector{Float64},h::Matrix{Float64},
        trace::Vector{LSTrace},i::Int,(mode::String,))

Return Newton's descend direction ``h^(-1)*g``.

# args
- `mode::String`: "cholesky" or (literally anything else, haven't implement more modes yet).
        Falls back to Julia's matrix inversion if ``h`` is not Cholesky decomposable.
"""
function solverNewton(booster::Booster,hist::Vector{State},freqs::Vector{Float64},
        objFunction::Callback,
        p::Vector{Float64},g::Vector{Float64},h::Matrix{Float64},
        trace::Vector{LSTrace},i::Int,(mode,)::Tuple{String,})
    
    if lowercase(mode) == "cholesky"
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

"""
    SolverNewotn(mode)

Callback for Newton's descend solver in [`linesearch`](@ref). mode is either "cholesky"
or not. See [`solverNewton`](@ref).
"""
SolverNewton(mode::String) = Callback(solverNewton,(mode,))


"""
    solverBFGS(booster::Booster,hist::Vector{State},freqs::Vector{Float64},
        objFunction::Callback,p,g,h,trace,i,(h0,))

Return quasi-Newtonian [`BFGS`](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher\
%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) descend direction.

# args
- `h0`: Hessian-like matrix for initialization of BFGS approximate Hessian.
"""
function solverBFGS(booster::Booster,hist::Vector{State},freqs::Vector{Float64},
        objFunction::Callback,p,g,h,trace,i,(h0,))
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

"""
    SolverBFGS(h0)

Callback for BFGS descend solver in [`linesearch`](@ref). h0 is Hessian for initialization.
See [`solverBFGS`](@ref).
"""
SolverBFGS(h0) = Callback(solverBFGS,(h0,))




"""
    solverHybrid(booster::Booster,hist::Vector{State},freqs::Vector{Float64},
        objFunction::Callback,p::Vector{Float64},g::Vector{Float64},h::Matrix{Float64},
        trace::Vector{Dragoon.LSTrace},i::Int,
        (mode,ϵls,αtest,ntest)::Tuple{String,Real,Real,Int})

Return Newton's descend direction if it actually descends, fallback to steepest descend if
not. Tests Newton step forwards and backwards similar to discrete linesearch.

# args
- `mode::String`: Newton solver mode, see [`solverNewton`](@ref).
- `ϵls::Real`: Minimum required descend for a single Newton test step.
- `αtest::Real`: Length of test step.
- `ntest::Int`: Amount of test steps to perform.
"""
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



"""
    SolverHybrid(mode::String,ϵls::Real,αtest::Real,ntest::Int)

Callback for Hybrid descend solver in [`linesearch`](@ref). See [`solverBFGS`](@ref).
"""
SolverHybrid(mode::String,ϵls::Real,αtest::Real,ntest::Int) =
    Callback(solverHybrid,(mode,ϵls,αtest,ntest))
