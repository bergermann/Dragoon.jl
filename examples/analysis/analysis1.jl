
using Pkg; Pkg.update()
using Dragoon, BoostFractor
using HDF5, DataFrames

include("tools.jl");

# path = getPath(1e-3,100_000,22.025e9,50e6,10,20,24.0,0.0);
# path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0,"NM2","2024_05_29-15_05_16");

# data = prepareDataAll1d(path,-10_000);
data = prepareDataAll1d(getPath(),-14_500);
# data = prepareData1d(path,-14_000);





initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);

sortData!(data)
b = best(data);
# showFields(data,data.freqs,freqsplot)
# showDist(data[2])

# showQuality(data)

# showDist(data,1000; xlabel="disc index", ylabel="d_i [mm]")
# hline!([initdist]*1e3)

# plot(data.freqs/1e9,data.boost; xlabel="frequency [GHz]",ylabel="boost factor β^2")

c = showClusters(data,50,200);
bc = [best(data[findall(isequal(i),c.assignments)]) for i in axes(c.centers,2)]
# scanSingleDiscs(best(data).dist,22.025e9,50e6,50,250,100)

showDistribution(data,d0)

datanm = data[findall(isequal(:nm),data.tags)]; sortData!(datanm)
datasa = data[findall(isequal(:sa),data.tags)]; sortData!(datasa)
datals = data[findall(isequal(:ls),data.tags)]; sortData!(datals)

showDistribution(datanm,d0)
showDistribution(datasa,d0)
showDistribution(datals,d0)



out = findOutliers(data,4.5e-3; showdistribution=true); length(out)
showDist(out)
showFields(out,data.freqs,freqsplot)



wiggle(data[1].pos[:,1],1e-6,10000,freqsplot,(22e9,22.05e9))
wiggle(data[2].pos[:,1],1e-6,10000,freqsplot,(22e9,22.05e9))

wigglewiggle(data[1].pos[:,1],collect(1:10)*1e-6,10000,freqsplot,(22e9,22.05e9))
wigglewiggle(data[2].pos[:,1],collect(1:10)*1e-6,10000,freqsplot,(22e9,22.05e9))


wiggle(data[1].pos[:,1],5e-6,10000,freqsplot,(22e9,22.05e9))
wiggle(data[2].pos[:,1],5e-6,10000,freqsplot,(22e9,22.05e9))

