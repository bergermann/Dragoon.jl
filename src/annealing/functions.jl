




function printSAIter(booster::Booster,obj::Float64,objsol::Float64,τ::Float64,iter::Int)
    if hasproperty(booster,:startingtime)
        println("Iter: ",iter,", timestamp: ",canonicalize(
            round(booster.timestamp-booster.startingtime,Second)))
    else
        println("Iter: ",iter,", timestamp: ",canonicalize(
            round(booster.timestamp-DateTime(0),Second)))
    end

    println("Iter finished. Objective value current:  ",round(obj; digits=3),"\n")
    println("               Objective value solution: ",round(objsol; digits=3),"\n")
    println("               Temperature:              ",round(τ; digits=3),"\n")
end





mutable struct SATrace
    x::Array{Float64}
    obj::Float64
    xsol::Array{Float64}
    objsol::Float64
    τ::Float64
    iter::Int
    t::DateTime
    T::Float64

    function SATrace(x,obj,xsol,objsol,τ,iter,t,T)
        new(x,obj,xsol,objsol,τ,iter,t,T)
    end

    function SATrace()
        new([0],0,0,0,DateTime(0),0)
    end

    function SATrace(x,obj,τ,iter)
        new(x,obj,τ,iter,DateTime(0),0)
    end    
end

function findNeighbour(booster::Booster,rmax::Float64)
    x_ = 2*rand(Float64,booster.ndisk) .- 1
    x_ /= pNorm(x_)

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


