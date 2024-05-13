using JLD2, Plots, ColorSchemes, Clustering, BoostFractor, Dragoon

include("tools.jl")

path = getPath(1e-6,100_000,22.025e9,50e6,10,20,24.0,0.0);

data = prepareDataAll1d(path,-12_000);

# data = prepareData1d(path,-12_000);

showQuality(data,-14000)

showDist(data,1000)
c = showClusters(data,50,100);

