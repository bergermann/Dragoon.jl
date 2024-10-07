
using Pkg; Pkg.update()


using Dragoon
using Plots, Plots.Measures
using Dates
using JLD2

include(joinpath(pwd(),"examples\\analysis\\tools\\tools.jl"));

@load "examples\\full_20_24.0_6.0e-5.jld2"

f = 25

freqs = genFreqs(f*1e9+25e6,100e6; n=50);
freqsplot = genFreqs(f*1e9+25e6,150e6; n=1000);

booster = AnalyticalBooster(P0[f]; tand=6e-5)
obj = ObjAnalytical
hist = initHist(booster,2*(booster.ndisk^2),freqs,obj);

plot(freqsplot/1e9,getBoost1d(booster,freqsplot))


ref0 = getRef1d(booster,freqs)
plot(freqs/1e9,abs.(ref0))


d0 = findpeak1d((freqs[1]+freqs[end])/2,booster.ndisk;
    eps=booster.epsilon,tand=booster.tand,granularity=10_000,deviation=0.3)
move(booster,dist2pos(ones(booster.ndisk)*d0); additive=false)

plot(freqsplot/1e9,getBoost1d(booster,freqsplot)/1e3)
plot(freqsplot/1e9,abs.(getRef1d(booster,freqsplot)))


move(booster,dist2pos(ones(booster.ndisk)*d0); additive=false)
nelderMead(booster,hist,freqs,
            1.,1+2/booster.ndisk,0.75-1/2booster.ndisk,1-1/booster.ndisk,1e-5,
            ObjRef1dS(ref0,x->x),
            InitSimplexRegular(1e-4),
            DefaultSimplexSampler,
            UnstuckNew(InitSimplexRegular(1e-4),true,5);
            maxiter=Int(5e1),
            showtrace=true,
            showevery=Int(1e0),
            unstuckisiter=true,);

            
plot(freqsplot/1e9,getBoost1d(booster,freqsplot))



move(booster,dist2pos(ones(booster.ndisk)*d0); additive=false)
trace = simulatedAnnealing(booster,hist,freqs,
            100e-6,
            ObjRef1dS(ref0,x->x),
            ObjAnalytical,
            UnstuckDont;
            maxiter=Int(10001),
            nreset=500,
            nresetterm=10,
            showtrace=true,
            showevery=1000,
            unstuckisiter=true,
            traceevery=100);
            
plot(freqsplot/1e9,getBoost1d(booster,freqsplot))