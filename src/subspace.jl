
using LinearAlgebra: eigvals, eigen, eigvecs


function subspace(booster,hist,freqs,lb,ub,n,derivator)
    @assert length(lb) == length(ub) == booster.ndisk "Bounds need same length as disks."

    # G = zeros(Float64,n,booster.ndisk)
    M = zeros(booster.ndisk,booster.ndisk)
    
    g = zeros(Float64,booster.ndisk); h = zeros(0,0)

    for i in 1:n
        move(b,rand(Float64,booster.ndisk)*(ub-lb)+lb)

        derivator.func(g,h,booster,hist,freqs,objFunction,derivator.args)

        M += a*a'
    end

    M ./= n

    E = eigen(M; sortby=x->-abs(x))

    return
end
