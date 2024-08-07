
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
        0.0069424515103718450,
        0.0070364724794998870,
        0.0078732484028712010,
        0.0070004748835278420,
        0.0073614019954116670,
        0.0071549044978148910,
        0.0070875765947395550,
        0.0074426096458853840,
        0.0073301995131451340,
        0.0074220155045692520,
        0.0071730390572375410,
        0.0070696702277422010,
        0.0073113174166905695,
        0.0065711615775509440,
        0.0095807997342139400,
        0.0048028172495510070,
        0.0079482239556944450,
        0.0064862237924530890,
        0.0076375922093132050,
        0.0068707736863909690,
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
                ObjRefSquare(ref0),
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
        "NM1REF"
    )

if !isdir(path)
    mkpath(path)
end

println("saving to $(joinpath(path,"$(date).jld2"))")

@save joinpath(path,"$(date).jld2") data sigx s seed T