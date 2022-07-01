export solverSteepestDescent, solverNewton, solverBFGS

import LinearAlgebra: cholesky

#solvers to determine directions

function solverSteepestDescent(p,g,h,trace,i)
    p[:] = -g
end

function solverNewton(p,g,h,trace,i,mode)
    if mode == "cholesky"
        C = cholesky(h)

        p[:] = inv(C.U)*inv(C.L)*g
    else
        p[:] = inv(h)*g
    end
end

function solverBFGS(p,g,h,trace,i,h0)
    if i == 1
        trace[i].h = h0

        p[:] = inv(h0)*g
    else
        y = g - trace[i-1].g
        s = trace[i].x - trace[i-1].x
        h_ = trace[i-1].h
        h = h_ + y*y'/(y'*s) - h_*s*s'*h_'/(s'*h_*s)
        trace[i].h = h

        p[:] = inv(h)*g
    end
end
