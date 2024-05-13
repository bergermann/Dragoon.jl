using JLD2, Plots, ColorSchemes, Clustering, BoostFractor, Dragoon

include("tools.jl")

path = joinpath(
    "examples",
    "analysis",
    "optimization data",
    "1.0e-6_100000_2.2025e10_5.0e7_10_20_24.0_0.0",
    "NM1",
    "2024_05_06-08_12_59.jld2"
);

data = prepareData1d(path,-12_000);

histogramB(data.obj)
histogramB(data.optdist)
histogramB(data.opttime)
histogramB(data.runtime)

showDist(data,1000)
showClusters(data,50,100);

showQuality(data,-14000)

plot(freqs,boosts[1:10_000,:]'; legend=false,label="")