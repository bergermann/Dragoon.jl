
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

const initdist = findpeak1d(s.f0,s.ndisk)
pos0 = dist2pos(ones(s.ndisk)*initdist);

booster = AnalyticalBooster(initdist; ndisk=s.ndisk,ϵ=s.eps,tand=s.tand)

@everywhere using Dragoon, Random 

@everywhere begin
    const sigx = $sigx
    const freqs = $freqs
    const pos0 = $pos0
    const booster = $booster

    const ref0 = [
        0.921515054155135100 + 0.38834263861365180im,
       -0.477337960337681830 + 0.87871979129906400im,
       -0.990553638519487400 + 0.13712581527853937im,
       -0.917166074070469500 - 0.39850519767524295im,
       -0.660353247163839800 - 0.75095511780676720im,
       -0.272694156192703760 - 0.96210077288107000im,
        0.228479273044231430 - 0.97354877730352020im,
        0.727006844461289600 - 0.68663021205481320im,
        0.990485300895356000 - 0.13761856237514040im,
        0.906491430910099100 + 0.42222421257735720im,
    ]
end

const pids = ParallelUtilities.workers_myhost()

data = SharedArray{Float64}((Nsig,s.ndisk+4); pids=pids)
T = SharedArray{Float64}(Nsig; pids=pids)

const seed = rand(UInt)

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
