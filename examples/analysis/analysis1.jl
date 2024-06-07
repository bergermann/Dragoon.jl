using JLD2, Plots, ColorSchemes, Clustering, BoostFractor, Dragoon, Statistics, Peaks

include("tools.jl")

# path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0);
path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0,"NM2","2024_05_29-15_05_16");
# path = joinpath("examples","analysis","optimization data","2024_05_27-15_21_00.jld2");

# data = prepareDataAll1d(path,-10_000);
data = prepareData1d(path,-14_000);

initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk);
sortData!(data)

# mean()

showQuality(data,0)

showDist(data,1000; xlabel="disc index", ylabel="d_i [mm]")
hline!([initdist]*1e3)

plot(data.freqs/1e9,data.boost; xlabel="frequency [GHz]",ylabel="boost factor β^2")

c = showClusters(data,100,200);

freqsplot = genFreqs(22.025e9,150e6; n=1000);


b = best(data);

plot(data.freqs/1e9,b.boost)
scatter(b.dist)

scanSingleDiscs(b.dist,data.s.f0,data.s.df,10,3_000,1_000)





disc = 15
l = λ(data.s.f0); n = 300_000
B = zeros(10,n)
d_ = copy(b.dist); d_[disc] = 0
d = zeros(20); d[disc] = 1

for i in 1:n
    B[:,i] = boost1d(d_+d*l*i/100_000,data.freqs)
end

obj = reshape(minimum(B,dims=1),(n,));

p = findmaxima(obj,40000)

plt = plot((1:n)/100_000,obj; xlabel="d_$disc/λ",ylabel="minimum boost")
vline!(p.indices/100_000)
display(plt)

pp = p.indices*l/100_000
display((pp[2:end] - pp[1:end-1])/(l/2))

plot(data.freqs/1e9,B[:,p.indices])