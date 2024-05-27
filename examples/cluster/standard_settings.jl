
mutable struct Settings
    f0::Float64
    df::Float64
    nf::Int
    ndisk::Int
    eps::Float64
    tand::Float64

    function Settings()
        new(
            22.025e9,   # f0::Float64
            50e6,       # df::Float64
            10,         # nf::Int
            20,         # ndisk::Int
            24.0,       # eps::Float64
            0.0,        # tand::Float64
        )
    end

    function Settings(f0,df,nf,ndisk,eps,tand)
        new(f0,df,nf,ndisk,eps,tand)
    end
end



function parseArgs(args)
    sigx = parse(Float64,args[1])
    Nsig = parse(Int,args[2])
    
    println("σ: ",sigx)
    println("N: ",Nsig,"\n")

    s = Settings()
    
    for (i,arg) in enumerate(args[3:min(8,length(args))])
        if arg == "_" || arg == "*"
            continue
        else
            setfield!(s,i,parse(fieldtype(Settings,i),arg))
        end
    end
    
    println(s)
    println("\n")

    return sigx, Nsig, s, nothing
end



import Base: println
function Base.println(s::Settings)
    println("--- Settings ---")
    println("Center frequency: $(s.f0)")
    println("Span frequency:   $(s.df)")
    println("Frequency points: $(s.nf)")
    println("Disc amount:      $(s.ndisk)")
    println("Disc epsilon:     $(s.eps)")
    println("Disc loss (tand): $(s.tand)")
end

function printOutput(data,T,ndisk)
    println()
    println("Best objective value:\n$(minimum(data[:,ndisk+1]))")
    println("Worst objective value:\n$(maximum(data[:,ndisk+1]))")
    println()
    println("Longest optimization time:\n$(maximum(T))")
    println("Shortest optimization time:\n$(minimum(T))")
    println("Average optimization time:\n$(sum(T)/length(T))")
    println()

    return
end