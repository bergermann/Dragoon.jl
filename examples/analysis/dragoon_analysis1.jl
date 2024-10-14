
using Pkg; Pkg.update()


using Dragoon
using Plots, Plots.Measures
using Dates
using JLD2

include(joinpath(pwd(),"examples\\analysis\\tools\\tools.jl"));

@load "examples\\full_20_24.0_6.0e-5.jld2"

f = 22

freqs = genFreqs(f*1e9+25e6,50e6; n=10);
freqsplot = genFreqs(f*1e9+25e6,150e6; n=1000);

booster = AnalyticalBooster(P0[f]; tand=6e-5)
obj = ObjAnalytical
hist = initHist(booster,2*(booster.ndisk^2),freqs,obj);

plot(freqsplot/1e9,getBoost1d(booster,freqsplot))
minimum(getBoost1d(booster,freqs))



move(booster,P0[22]; additive=false);
O1,P1,F1,s1 = dragoon(booster,hist,50e6,5e6,
        obj,UnstuckDont;
        fmin=10.0e9,fmax=22.05e9,nfreqs=10,niters=2000,
        scalerange=(0.5,1.5),scalesteps=1000,
        preoptimize=false,reverse=true)
f1 = [(f[1]+f[end])/2 for f in F1];



p3 = plot(collect(10:1:100).+0.025,-B0/1e3,label="scratch",seriestype=:scatter,
    xlabel="Frequency [GHz]",ylabel="Objective Value × 10³",markersize=2)
plot!(f1/1e9,-O1/1e3; label="rescaling 1",c=:blue,lw=2)
plot(f1/1e9,-O1/1e3; label="rescaling 1",c=:blue,lw=2)


plot(s1)

p4 = plot(collect(10:1:100).+0.025,-B0/1e3,label="scratch",seriestype=:scatter,
    xlabel="Frequency [GHz]",ylabel="Objective Value × 10³",markersize=2,legend=false)
for f_ in 45:50
    println(f_)
    move(booster,P0[f_]; additive=false);
    O1,P1,F1,s1 = dragoon(booster,hist,50e6,5e6,
        obj,UnstuckDont;
        fmin=f_*1e9,fmax=(f_+1)*1e9,nfreqs=10,niters=3000,
        scalerange=(1.1,1.3),scalesteps=1000,
        preoptimize=false,reverse=false)
    f1 = [(f[1]+f[end])/2 for f in F1];
    plot!(p4,f1/1e9,-O1/1e3; c=:blue,lw=2,label="")
end
p4


# savefig(p5,"rescale.svg")

plotRescale(booster,P0[22],22.025e9,50e6,200e6,1.175,-5,5)


p6 = plot(; legend=false); p7 = plot(; legend=false)

# for i in cat(10:27,33:57; dims=1)
for i in 33:57
    f = genFreqs(i*1e9+25e6,150e6; n=100)
    move(booster,P0[i])
    # ref = abs.(getRef1d(booster,f))
    boost = normalize_range(getBoost1d(booster,f))
    ref = normalize_range(abs.(getRef1d(booster,f)))

    # plot!(p6,log.(boost))
    # plot!(p7,log.(ref))
    plot!(p6,boost)
    plot!(p7,ref)
end

display(p6); display(p7)

