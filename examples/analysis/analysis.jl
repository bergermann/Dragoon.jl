using JLD2, Plots, ColorSchemes, Clustering, BoostFractor, Dragoon

include("tools.jl")

path = getPath(1e-5,1_000,22.025e9,50e6,10,20,24.0,0.0,"NM1ref");

data = prepareDataAll1d(path,0.1);

# data = prepareData1d(path,-12_000);

showQuality(data,0.1)

showDist(data,1000)
c = showClusters(data,10,100);

