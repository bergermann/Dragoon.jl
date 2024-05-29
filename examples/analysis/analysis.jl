using JLD2, Plots, ColorSchemes, Clustering, BoostFractor, Dragoon, Statistics

include("tools.jl")

path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0);
# path = joinpath("examples","analysis","optimization data","2024_05_27-15_21_00.jld2");

data = prepareDataAll1d(path,-10_000);
# data = prepareData1d(path,0);

initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk);
sortData!(data)

mean()

# data = prepareDataAll1d(path,-12_000);

# showQuality(data,0)

showDist(data,1000)
hline!([initdist])

c = showClusters(data,100,200);

freqsplot = genFreqs(22.025e9,150e6; n=1000);

plot(data.freqs/1e9,)

b = best(data);

plot(data.freqs/1e9,b.boost)
scatter(b.dist)