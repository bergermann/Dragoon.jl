using JLD2, Plots, ColorSchemes, Clustering, BoostFractor, Dragoon

include("tools.jl")

# path = getPath(1e-5,1000,22.025e9,50e6,50,20,24.0,0.0);
path = joinpath("examples","analysis","optimization data","2024_05_27-15_21_00.jld2");

# data = prepareDataAll1d(path);
data = prepareData1d(path,0);

# data = prepareDataAll1d(path,-12_000);

showQuality(data,0)

showDist(data,1000)
c = showClusters(data,7,200);

freqsplot = genFreqs(22.025e9,150e6; n=1000);

plot(data.freqs/1e9,)

b = best(data)

plot(data.freqs/1e9,b.boost)
scatter(b.dist)