
using Distributed, ParallelUtilities, SharedArrays, JLD2, Dragoon
@everywhere using Dragoon, Random

println("Available processors: ",nprocs())
println("Available workers:    ",nworkers(),"\n")

include("standard_settings.jl")



function main(args)
    sigx, Nsig, s, _ = parseArgs(args)

    freqs = genFreqs(s.f0,s.df; n=s.nf)

    initdist = findpeak1d(s.f0,s.ndisk)
    pos0 = dist2pos(ones(s.ndisk)*initdist);

    booster = AnalyticalBooster(initdist; ndisk=s.ndisk,ϵ=s.eps,tand=s.tand)
    booster.wavelength = λ(s.f0)
    booster.mindist = 1e-3
    
    dist0 = [
        0.0069864996648397000,
        0.0072560130634651010,
        0.0072587766881647050,
        0.0072180129444167985,
        0.0071810463115898690,
        0.0072643048247736110,
        0.0072005219024528160,
        0.0072281469176733665,
        0.0072039644723843060,
        0.0072056913829301470,
        0.0072293906115958880,
        0.0071775904582249700,
        0.0072014504364220360,
        0.0073195415726201870,
        0.0071448147947109080,
        0.0073325781037591470,
        0.0071798833031334620,
        0.0073681846490424320,
        0.0066949319976553630,
        0.0074786250621039780,
    ]

    ref0 = ref1d(dist0,freqs; eps=s.eps,tand=s.tand)

    @everywhere begin
        sigx = $sigx
        freqs = $freqs
        pos0 = $pos0
        booster = $booster

        ref0 = $ref0
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
            
            trace, term = nelderMead(booster,hist,freqs,
                1.,1+2/booster.ndisk,0.75-1/(2*booster.ndisk),1-1/(booster.ndisk),1e-6,
                ObjRef1dSquare(ref0),
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

    return data, sigx, Nsig, s, seed, T
end

data, sigx, Nsig, s, seed, T = main(ARGS)

date = getDateString()
path = joinpath(
        "optimization data",
        "$(sigx)_$(Nsig)_$(s.f0)_$(s.df)_$(s.nf)_$(s.ndisk)_$(s.eps)_$(s.tand)",
        "NM2REF"
    )

if !isdir(path)
    mkpath(path)
end

println("saving to $(joinpath(path,"$(date).jld2"))")

@save joinpath(path,"$(date).jld2") data sigx s seed T