
using Dragoon, BoostFractor
using JLD2, HDF5, DataFrames

include("tools.jl");

# path = getPath(1e-3,100_000,22.025e9,50e6,10,20,24.0,0.0);
path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0,"NM2","2024_05_29-15_05_16");

# data = prepareDataAll1d(path,-10_000);
# data = prepareDataAll1d(getPath(),-14_000);
data = prepareData1d(path,-14_000);



function saveData(data)
    if is
        df = DataFrame()

        for i in 1:data.s.ndisk
            df[!,"d$i"] = data.dist[i,:]
        end
        df[!,"obj"] = data.obj
    end

    sort!(df,[:obj])
    unique!(df)

    h5write("test.h5","data",df)

    return
end













initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=200);

sortData!(data)
b = best(data)

# showQuality(data)

# showDist(data,1000; xlabel="disc index", ylabel="d_i [mm]")
# hline!([initdist]*1e3)

# # plot(data.freqs/1e9,data.boost; xlabel="frequency [GHz]",ylabel="boost factor β^2")

# # c = showClusters(data,100,200);
# # scanSingleDiscs(best(data).dist,22.025e9,50e6,50,250,100)

# showDistribution(data)


showFields(b.pos,freqsplot)