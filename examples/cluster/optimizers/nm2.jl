
using Distributed, ParallelUtilities, SharedArrays, JLD2, Dragoon

println("Available processors: ",nprocs())
println("Available workers:    ",nworkers(),"\n")

include("standard_settings.jl")



function main(args)
    sigx, Nsig, s, _ = parseArgs(args)

    freqs = genFreqs(s.f0,s.df; n=s.nf)

    booster = AnalyticalBooster(1e-3; ndisk=s.ndisk,ϵ=s.eps,tand=s.tand)

    @everywhere using Dragoon, Random 

    @everywhere begin
        sigx = $sigx
        freqs = $freqs
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

            λ = 299792458.0/s.f0

            p0 = dist2pos(randn(booster.ndisk)*sigx*λ)

            move(booster,p0; additive=false)

            hist = initHist(booster,100,freqs,ObjAnalytical)
            booster.summeddistance = 0.
            
            trace, term = nelderMead(booster,hist,freqs,
                1.,1+2/booster.ndisk,0.75-1/(2*booster.ndisk),1-1/(booster.ndisk),1e-6,
                ObjAnalytical,
                InitSimplexRegular(5e-5),
                DefaultSimplexSampler,
                UnstuckNew(InitSimplexRegular(5e-5),true,-10000);
                maxiter=2000,
                traceevery=typemax(Int),
                showtrace=false,
                unstuckisiter=true,
                resettimer=true,
                returntimes=true)

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
            "rand_$(Nsig)_$(s.f0)_$(s.df)_$(s.nf)_$(s.ndisk)_$(s.eps)_$(s.tand)",
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