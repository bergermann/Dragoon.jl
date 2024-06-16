using Plots, ColorSchemes, Clustering, BoostFractor, Dragoon, Statistics, Peaks
using JLD2, HDF5, DataFrames

include("tools.jl");

path = getPath(1e-3,100_000,22.025e9,50e6,10,20,24.0,0.0);
# path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0,"NM2","2024_05_29-15_05_16");

# data = prepareDataAll1d(path,-10_000);
data = prepareDataAll1d(getPath(),-14_000);
# data = prepareData1d(path,-14_000);

initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);

sortData!(data)
b = best(data)

showQuality(data)

showDist(data,1000; xlabel="disc index", ylabel="d_i [mm]")
hline!([initdist]*1e3)

# plot(data.freqs/1e9,data.boost; xlabel="frequency [GHz]",ylabel="boost factor β^2")

# c = showClusters(data,100,200);

showDistribution(data)


function wiggle(pos::Vector{Float64},sigx::Real,n::Int,
                freqsplot::Vector{Float64},bounds::Tuple{Float64,Float64})

    d = pos2dist(pos)
    b0 = boost1d(d,freqsplot)
    obj0 = minimum(b0[@.(bounds[1] < freqsplot < bounds[2])])

    db = zeros(Float64,length(freqsplot),n)
    obj = zeros(Float64,n)

    for i in 1:n
        db[:,i] = boost1d(d+sigx*randn(length(pos)),freqsplot)
        obj[i] = minimum(db[@.(bounds[1] < freqsplot < bounds[2]),i])
    end

    db .-= b0

    q = collect(0:0.1:1)
    Q = zeros(length(q),length(freqsplot))
    for i in axes(Q,2)
        Q[:,i] = quantile(db[i,:],q)
    end

    p1 = plot(freqsplot/1e9,b0; xlabel="Frequency [GHz]",label="original")
    plot!(p1,freqsplot/1e9,b0+Q[div(length(q)-1,2),:],label="50% quant")

    for i in 1:div(length(q)-1,2)
        plot!(p1,freqsplot/1e9,b0+Q[i,:]; linewidth=0,
            fillrange=b0+Q[end-i+1,:],fillalpha=0.3,fillcolor=1,label=""
            )
    end

    p2 = histogram(obj)
    vline!(p2,[obj0])

    display(p1)
    display(p2)

    return Q
end




function saveData(data)
    df = DataFrame()

    for i in 1:data.s.ndisk
        df[!,"d$i"] = data.dist[i,:]
    end
    df[!,"obj"] = data.obj

    sort!(df,[:obj])
    unique!(df)
end


h5write("test.h5","data",df)