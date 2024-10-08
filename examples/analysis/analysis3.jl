
using Pkg; Pkg.update()
using Dragoon, BoostFractor
using HDF5, DataFrames

include("tools.jl");

data = prepareDataAll1d(getPath(),-15_000);

initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);

sortData!(data)
b = best(data);

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



wiggle(data[1].pos[:],1e-6,100,freqsplot,(22e9,22.05e9))
wiggle(data[2].pos[:],1e-6,100,freqsplot,(22e9,22.05e9))

wigglewiggle(data[1].pos[:],collect(1:10)*1e-6,10000,freqsplot,(22e9,22.05e9))
wigglewiggle(data[2].pos[:],collect(1:10)*1e-6,10000,freqsplot,(22e9,22.05e9))





r1 = ref1d(data[1].dist[:],freqsplot; tand=1e-4);
r2 = ref1d(data[2].dist[:],freqsplot; tand=1e-4);

b1_ = boost1d(data[1].dist[:],freqsplot; tand=0.);
b2_ = boost1d(data[2].dist[:],freqsplot; tand=0.);
b1 = boost1d(data[1].dist[:],freqsplot; tand=1e-5);
b2 = boost1d(data[2].dist[:],freqsplot; tand=1e-5);

plot(freqsplot/1e9,b1)
plot!(freqsplot/1e9,b2)
plot!(freqsplot/1e9,b1_)
plot!(freqsplot/1e9,b2_)

plot(freqsplot/1e9,abs.(r1))
plot!(freqsplot/1e9,abs.(r2))

plot(data.freqs/1e9,real.(data[1].ref))
plot!(data.freqs/1e9,imag.(data[1].ref))

plot(data.freqs/1e9,real.(data[2].ref))
plot!(data.freqs/1e9,imag.(data[2].ref))