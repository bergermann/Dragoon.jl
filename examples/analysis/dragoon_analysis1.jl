
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
        fmin=10.0e9,fmax=22.05e9,nfreqs=10,niters=5000,
        scalerange=(0.5,1.5),scalesteps=2000,
        preoptimize=false,reverse=true)
f1 = [(f[1]+f[end])/2 for f in F1];



p3 = plot(collect(10:1:100).+0.025,-B0/1e3,label="scratch",seriestype=:scatter,
    xlabel="Frequency [GHz]",ylabel="Objective Value × 10³",markersize=3)
plot!(f1/1e9,-O1/1e3; label="rescaling 1",c=:blue,lw=2)
# plot(f1/1e9,-O1/1e3; label="rescaling 1",c=:blue,lw=2)


plot(s1)



p4 = plot(collect(10:1:100).+0.025,-B0/1e3,label="scratch",seriestype=:scatter,
    xlabel="Frequency [GHz]",ylabel="Objective Value × 10³",markersize=3,legend=true)
for f_ in 10:100
    println(f_)
    move(booster,P0[f_]; additive=false);
    O1,P1,F1,s1 = dragoon(booster,hist,50e6,5e6,
        obj,UnstuckDont;
        fmin=f_*1e9,fmax=(f_+1)*1e9,nfreqs=10,niters=1000,
        scalerange=(1.0,2.0),scalesteps=2000,
        preoptimize=false,reverse=false)
    f1 = [(f[1]+f[end])/2 for f in F1];
    plot!(p4,f1/1e9,-O1/1e3; c=:blue,lw=2,label="",alpha=0.5)
end
plot(p4,[],label="rescaling",c=:blue,lw=2,alpha=0.5,legend=true)


# savefig(p5,"rescale.svg")

# plotRescale(booster,P0[22],22.025e9,50e6,200e6,1.175,-5,5)


p6 = plot(; xlabel="Frequency Index i",ylabel="Unnormalised Boost",legend=false,title="Area 1");
p7 = plot(; xlabel="Frequency Index i",ylabel=L"Unnormalised Reflectivity $|R|$",legend=false,title="Area 1");

# for i in cat(10:27,33:57; dims=1)
for i in 10:27
    f = genFreqs(i*1e9+25e6,150e6; n=100)
    move(booster,P0[i])
    boost = getBoost1d(booster,f)
    ref = abs.(getRef1d(booster,f))
    # boost = normalize_range(getBoost1d(booster,f))
    # ref = normalize_range(abs.(getRef1d(booster,f)))

    # plot!(p6,log.(boost))
    # plot!(p7,log.(ref))
    plot!(p6,boost)
    plot!(p7,ref)
end

display(p6); display(p7)

plot!(p3,collect(10:27),-B0[1:18]/1e3; label="area 1",lw=3,c=:red,alpha=0.5)
plot!(p3,collect(33:56),-B0[24:47]/1e3; label="area 2",lw=3,c=:green,alpha=0.5)
plot!(p3,collect(67:80),-B0[58:71]/1e3; label="area 3",lw=3,c=:blue,alpha=0.5)
plot!(p3,collect(95:100),-B0[86:end]/1e3; label="area 3",lw=3,c=:yellow,alpha=0.5)


