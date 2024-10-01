
using Pkg; Pkg.update()


using Dragoon
using Plots, Plots.Measures
using Dates
using JLD2

include(joinpath(pwd(),"analysis\\tools\\tools.jl"));

@load "full_20_24.0.jld2"

n = 20; f = 22

freqs = genFreqs(f*1e9+25e6,50e6; n=10)
freqsplot = genFreqs(f*1e9+25e6,150e6; n=1000)

booster = AnalyticalBooster(P0[f]; tand=6e-5)

obj = ObjAnalytical

hist = initHist(booster,2*(booster.ndisk^2),freqs,obj);

plot(freqsplot/1e9,getBoost1d(booster,freqsplot))
minimum(getBoost1d(booster,freqs))



move(booster,P0[10]; additive=false)

O1,P1,F1,s1,s_1 = dragoon(booster,hist,50e6,5e6,

        obj,UnstuckDont;
        fmin=10e9,fmax=20.05e9,nfreqs=10,
        scalerange=(1.0,1.3),scalesteps=1000,
        preoptimize=false,reverse=false)

f1 = [(f[1]+f[end])/2 for f in F1];

p3 = plot(collect(10:1:100).+0.025,-O0/1e3,label="scratch",seriestype=:scatter,
    xlabel="Frequency [GHz]",ylabel="Objective Value × 10³",markersize=2)


plot!(p3,f1/1e9,-O1/1e3; label="rescaling 1",c=:blue,lw=2)

plot(s_1)