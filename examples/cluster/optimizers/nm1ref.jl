
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

freqs = genFreqs(s.f0,s.df; n=s.nf)

initdist = findpeak1d(s.f0,s.ndisk)
pos0 = dist2pos(ones(s.ndisk)*initdist);

booster = AnalyticalBooster(initdist; ndisk=s.ndisk,ϵ=s.eps,tand=s.tand)

dist0 = [
    0.0072448589080291290,
    0.0070711997431429300,
    0.0071246144242139010,
    0.0073342065660182530,
    0.0072327454408213510,
    0.0072511235941711350,
    0.0071632103241744970,
    0.0072023665514860190,
    0.0073048489540810120,
    0.0071862085372421955,
    0.0070924346714088270,
    0.0073231298843525620,
    0.0072555088582221690,
    0.0072399407925294880,
    0.0072071841247186580,
    0.0072590448084468580,
    0.0072497272553200550,
    0.0070498810897688120,
    0.0069997371085350490,
    0.0071354746829380430,
]

ref0 = ref1d(dist0,freqs; eps=s.eps,tand=s.tand)

@everywhere using Dragoon, Random 

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

date = getDateString()
path = joinpath(
        "optimization data",
        "$(sigx)_$(Nsig)_$(s.f0)_$(s.df)_$(s.nf)_$(s.ndisk)_$(s.eps)_$(s.tand)",
        "NM1"
    )

if !isdir(path)
    mkpath(path)
end

@save joinpath(path,"$(date).jld2") data sigx s seed T
