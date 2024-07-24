
using Pkg; Pkg.update()
using Dragoon, BoostFractor
using HDF5, DataFrames

include("tools.jl");

# path = getPath(1e-3,100_000,22.025e9,50e6,10,20,24.0,0.0);
# path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0,"NM2","2024_05_29-15_05_16");

# data = prepareDataAll1d(path,-10_000);
data = prepareDataAll1d(getPath(),-14_000);
# data = prepareData1d(path,-14_000);





initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);

sortData!(data)
b = best(data);

# showQuality(data)

# showDist(data,1000; xlabel="disc index", ylabel="d_i [mm]")
# hline!([initdist]*1e3)

# plot(data.freqs/1e9,data.boost; xlabel="frequency [GHz]",ylabel="boost factor β^2")

c = showClusters(data,50,200);
bc = [best(data[findall(isequal(i),c.assignments)]) for i in axes(c.centers,2)]
# scanSingleDiscs(best(data).dist,22.025e9,50e6,50,250,100)

# showDistribution(data)


# showFields(b.pos,freqsplot)
# showFields(b.pos,[data.freqs[2],sum(data.freqs[[2,8]])/2,data.freqs[8]])


# display(showFields(bc[i].pos,data.freqs[2])[1])



# P = [showFields(bc[i].pos,data.freqs[2]) for i in eachindex(bc)]

# ylim1 = maximum([ylims(p[1][1])[2] for p in P])
# ylim2 = maximum([ylims(p[1][2])[2] for p in P])

# for i in eachindex(P)
#     ylims!(P[i][1][1],(-ylim1,ylim1))
#     ylims!(P[i][1][2],(-ylim2,ylim2))
# end

# for p in P; display(p[1]); end


out = findOutliers(data,7e-3; showdistribution=true)

# f = data.freqs[5]
f = 22.0e9+1e6

showFields(out,data.freqs,freqsplot)
