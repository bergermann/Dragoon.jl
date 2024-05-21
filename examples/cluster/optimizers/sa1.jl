
using Distributed, ParallelUtilities, SharedArrays, JLD2, Dragoon

println("Available processors: ",nprocs())
println("Available workers:    ",nworkers(),"\n")

include("standard_settings.jl")

sigx = parse(Float64,ARGS[1])
Nsig = parse(Int,ARGS[2])

println("σ: ",sigx)
println("N: ",Nsig,"\n")

@assert length(ARGS[3:end]) <= fieldcount(Settings) 

for (i,arg) in enumerate(ARGS[3:end])
    if arg == "_" || arg == "*"
        continue
    else
        setfield!(s,i,parse(fieldtype(Settings,i),arg))
    end
end

println(s)
println("\n")

function main(args)
    freqs = genFreqs(s.f0,s.df; n=s.nf)

    initdist = findpeak1d(s.f0,s.ndisk)
    pos0 = dist2pos(ones(s.ndisk)*initdist);

    booster = AnalyticalBooster(initdist; ndisk=s.ndisk,ϵ=s.eps,tand=s.tand)

    @everywhere using Dragoon, Random 

    @everywhere begin
        sigx = $sigx
        freqs = $freqs
        pos0 = $pos0
        booster = $booster
    end

    pids = ParallelUtilities.workers_myhost()

    data = SharedArray{Float64}((Nsig,s.ndisk+4); pids=pids)
    T = SharedArray{Float64}(Nsig; pids=pids)

    seed = rand(UInt)

    println("Total loop time:")
    @time out = @distributed (+) for i in collect(1:Nsig)
        t = @elapsed begin
            Random.seed!(seed+i)

            move(booster,pos0+randn(booster.ndisk)*sigx; additive=false)

            hist = initHist(booster,100,freqs,ObjAnalytical)
            booster.summeddistance = 0.
            
            trace, term = simulatedAnnealing(booster,hist,freqs,
                100e-6,
                TempLinear(100,1001),
                ObjAnalytical,
                UnstuckDont;
                maxiter=Int(1e5),
                nreset=500,
                nresetterm=10,
                showtrace=false,
                unstuckisiter=true,
                resettimer=true,
                traceevery=typemax(Int),
                returntimes=true);

            data[i,1:booster.ndisk] .= booster.pos
            data[i,booster.ndisk+1] = term[1]
            data[i,booster.ndisk+2] = term[2][1].value
            data[i,booster.ndisk+3] = term[2][2]
            data[i,booster.ndisk+4] = term[2][3].value
        end

        T[i] = t

        0
    end

    printOutput(data,T,booster.ndisk)

    date = getDateString()
    path = joinpath(
            "optimization data",
            "$(sigx)_$(Nsig)_$(s.f0)_$(s.df)_$(s.nf)_$(s.ndisk)_$(s.eps)_$(s.tand)",
            uppercase(@__FILE__)[1:end-3]
        )

    if !isdir(path)
        mkpath(path)
    end

    println("saving to $(joinpath(path,"$(date).jld2"))")

    @save joinpath(path,"$(date).jld2") data sigx s seed T

    return
end

main(ARGS)