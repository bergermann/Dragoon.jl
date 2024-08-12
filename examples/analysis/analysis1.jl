
using Pkg; Pkg.update()
using Dragoon, BoostFractor
using HDF5, DataFrames

include("tools.jl");

data = prepareDataAll1d(getPath(),-14_000);
sortData!(data)

initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);

b = best(data);
# histogram(-data.obj; xlabel="Minimum Boostfactor β²",title="Distribution of Converged States (20 Discs, ε=24)",
#     label=false,xlims=(floor(minimum(-data.obj)),ceil(maximum(-data.obj))))
histogram(data.obj; xlabel="Minimum Boostfactor β²",title="Distribution of Converged States (20 Discs, ε=24)",
    label=false,xlims=(floor(minimum(data.obj)),ceil(maximum(data.obj))))
# showFields(data,data.freqs,freqsplot)
# showDist(data[2])

# showQuality(data)

showDist(data,1000; xlabel="Disc Index", ylabel="d_i [mm]")
# hline!([initdist]*1e3)

# plot(data.freqs/1e9,data.boost; xlabel="Frequency [GHz]",ylabel="Boostfactor β²")

c = showClusters(data,50,200);
bc = [best(data[findall(isequal(i),c.assignments)]) for i in axes(c.centers,2)]
# scanSingleDiscs(best(data).dist,22.025e9,50e6,50,250,100)

showDistribution(data,d0)

datanm = data[findall(isequal(:nm),data.tags)]; sortData!(datanm)
datasa = data[findall(isequal(:sa),data.tags)]; sortData!(datasa)
datals = data[findall(isequal(:ls),data.tags)]; sortData!(datals)

showDist(datanm)
showDistribution(datanm,d0)
showDist(datasa)
showDistribution(datasa,d0)
showDist(datals)
showDistribution(datals,d0)



out = findOutliers(data,6e-3; showdistribution=true); length(out)
showDist(out)
showFields(out,data.freqs,freqsplot)

showFields(data[2],data.freqs,freqsplot)

freqsplot1 = genFreqs(22.025e9,150e6; n=201);
wiggle(data[1].pos[:,1],1e-6,100,freqsplot,(22e9,22.05e9));
wiggle(data[2].pos[:,1],1e-6,100,freqsplot,(22e9,22.05e9));

wigglewiggle(data[1].pos[:,1],collect(1:10)*1e-6,100,freqsplot,(22e9,22.05e9))
wigglewiggle(data[2].pos[:,1],collect(1:10)*1e-6,100,freqsplot,(22e9,22.05e9))

wigglewiggle(data[2].pos[:,1],collect(1:10)*1e-6,100,freqsplot,(22e9,22.05e9))

wiggle(data[1].pos[:,1],5e-6,100,freqsplot,(22e9,22.05e9));
wiggle(data[2].pos[:,1],5e-6,100,freqsplot,(22e9,22.05e9));

