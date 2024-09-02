
using Pkg; Pkg.update()
using Dragoon, BoostFractor
using HDF5, DataFrames

include("tools/tools.jl");

data = prepareDataAll1d(getPath(),-14_000);
sortData!(data)

initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);

b = best(data);
histogram(-data.obj; xlabel="Minimum Boostfactor β²",title="Distribution of Converged States",
    label=false,xlims=(floor(minimum(-data.obj)),ceil(maximum(-data.obj))))
    
showQuality(data)

showDist(data,100; ndiv=10)

showDistribution(data,d0)

datanm = data[findall(isequal(:nm),data.tags)]; sortData!(datanm)
datasa = data[findall(isequal(:sa),data.tags)]; sortData!(datasa)
datals = data[findall(isequal(:ls),data.tags)]; sortData!(datals)

showDist(datanm,100)
showDistribution(datanm,d0)
showDist(datasa,100)
showDistribution(datasa,d0)
showDist(datals,100)
showDistribution(datals,d0)


idxs = [2,100,1_000,10_000,100_000,300_000,690_000];

showFields(data[idxs],data.freqs,freqsplot)



out = findOutliers(data,6e-3; showdistribution=true); length(out)
ins = data[findall(i->!(data.obj[i] in out.obj),eachindex(data))]
showDist(ins)
showDist(out)
showFields(out,data.freqs,freqsplot)

showDistribution(out,d0)


freqsplot1 = genFreqs(22.025e9,150e6; n=201);
wiggle(data.pos[:,1],1e-6,10_000,freqsplot1,(22e9,22.05e9));
wiggle(data.pos[:,2],1e-6,10_000,freqsplot1,(22e9,22.05e9));

wiggle(data[1].pos[:,1],5e-6,10_000,freqsplot1,(22e9,22.05e9));
wiggle(data[2].pos[:,1],5e-6,10_000,freqsplot1,(22e9,22.05e9));

wigglewiggle(data.pos[:,1],collect(1:10)*1e-6,1_000,freqsplot1,(22e9,22.05e9))
wigglewiggle(data.pos[:,2],collect(1:10)*1e-6,1_000,freqsplot1,(22e9,22.05e9))


for i in eachindex(out)
    wiggle(out.pos[:,i],1e-6,10_000,freqsplot1,(22e9,22.05e9));
    wigglewiggle(out.pos[:,i],collect(1:10)*1e-6,1_000,freqsplot1,(22e9,22.05e9))
end

for i in eachindex(idxs)
    wiggle(data.pos[:,idxs[i]],1e-6,10_000,freqsplot1,(22e9,22.05e9));
    wigglewiggle(data.pos[:,idxs[i]],collect(1:10)*1e-6,1_000,freqsplot1,(22e9,22.05e9))
end


longs = data[findall(i->data.optdist[i]>5,eachindex(data))]; length(longs)
# histogram(longs.optdist; xlabel="Travel Distance [m]")
histogram(data.optdist; xlabel="Travel Distance [m]",legend=false,ylim=[0,10_000],bins=200); vline!([5,])
showDistribution(longs,d0)