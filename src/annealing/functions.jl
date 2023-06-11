


mutable struct SATrace
    x::Array{Float64}
    obj::Float64
    τ::Float64
    iter::Int
    t::DateTime
    T::Float64

    function LSTrace(x,obj,τ,iter,t,T)
        new(x,obj,τ,iter,t,T)
    end

    function LSTrace()
        new([0],0,0,0,DateTime(0),0)
    end

    function LSTrace(x,obj,τ,iter)
        new(x,obj,τ,iter,DateTime(0),0)
    end    
end

function findNeighbour(booster::Booster,rmax::Float64)
    x_ = 2*rand(Float64,booster.ndisk) .- 1
    x_ /= pNorm(x)

    return rmax*rand()*x_
end

function thermal(objx::Float64,objy::Float64,T::Float64)
    return exp((objx-objy)/T)
end

function thermal(booster::Booster,obj0::Float64,pos2::Vector{Float64},
        objFunction::Callback)

    move()
    obj0
    
    return
end


