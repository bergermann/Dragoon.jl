
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

data = SharedArray{Float64}((Nsig,s.ndisk+4); pids=ParallelUtilities.workers_myhost())

seed = rand(UInt)

@time out = @distributed (+) for i in collect(1:Nsig)

    Random.seed!(seed+i)

    move(booster,pos0+randn(booster.ndisk)*sigx; additive=false)

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

    0
end

date = getDateString()
path = joinpath(
        "optimization data",
        "$(sigx)_$(Nsig)_$(s.f0)_$(s.df)_$(s.nf)_$(s.ndisk)_$(s.eps)_$(s.tand)",
        "NM1"
    )

if !isdir(path)
    mkpath(path)
end

@save joinpath(path,"$(date).jld2") data sigx s seed