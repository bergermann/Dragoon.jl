using JLD2, Plots, ColorSchemes, Clustering, BoostFractor, Dragoon

include("tools.jl")

path = joinpath(
    "examples",
    "analysis",
    "data",
    "1.0e-6_100000_2.2025e10_5.0e7_10_20_24.0_0.0",
    "NM1",
    "2024_05_06-08_12_59.jld2"
);

@load String(path) data sigx s seed T

data = data[(data[:,s.ndisk+1] .<= -12000),:];

datad = zeros(Float64,size(data,1),s.ndisk);
for i in axes(data,1)
    datad[i,:] = pos2dist(data[i,1:s.ndisk])
end

histogram(data[:,21])
histogram(T)
maximum(T)

histogramB(data[:,21])


scatter(1:20,datad[1:10_000,:]'; legend=false,c=:blue,markersize=2)
plot(1:20,datad[1:10_000,:]'; legend=false,c=:blue,markersize=2)

k = 50
m = kmeans(datad[1:10000,:]',k)

colors = palette([:red,:blue],k);
scatter(1:20,datad[1:100,:]'; c=colors[m.assignments'],legend=false)

ylim = [0.95*minimum(datad),1.05*maximum(datad)]
for i in 1:k
    p1 = scatter(1:s.ndisk,datad[1:10000,:][(m.assignments .== i),:]';
        c=colors[i],legend=false,ylims=ylim)
    display(p1)

    p2 = scatter(freqsm)
end

boosts = zeros(Float64,size(datad,1),s.nf);
refs = zeros(ComplexF64,size(datad,1),s.nf);
freqs = genFreqs(s.f0,s.df; n=s.nf);

for i in axes(datad,1)
    boosts[i,:] = boost1d(datad[i,:],freqs; eps=s.eps,tand=s.tand)
    refs[i,:] = ref1d(datad[i,:],freqs; eps=s.eps,tand=s.tand)
end

plot(freqs,boosts[1:10_000,:]'; legend=false,label="")