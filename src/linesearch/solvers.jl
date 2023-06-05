### solvers to determine directions

export solverSteepestDescent, solverNewton, solverBFGS
export SolverSteep, SolverNewton, SolverBFGS

import LinearAlgebra: cholesky


# args = ()
function solverSteepestDescent(p,g,h,trace,i,args)
    p[:] = -g
end

const SolverSteep = Callback(solverSteepestDescent)



# args = (mode,)
function solverNewton(p,g,h,trace,i,args)
    if args[1] == "cholesky"
        C = cholesky(h)

        p[:] = inv(C.U)*inv(C.L)*g
    else
        p[:] = inv(h)*g
    end
end

const SolverNewton(mode) = Callback(solverNewton,(mode,))



# args = (h0,)
function solverBFGS(p,g,h,trace,i,args)
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
