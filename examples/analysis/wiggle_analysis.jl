
using Pkg; Pkg.update()

using Dragoon
using JLD2, Peaks

include("tools/tools.jl");
@load "examples\\full_20_24.0_6.0e-5.jld2"


D = Dict{Int,Data}()

for i in 10:100
    try
        threshold = B0[i-9]*0.8
        data = prepareDataAll1d(getPath(),threshold; f0=i*1e9+25e6,df=50e6,tand=6e-5);
        sortData!(data)

        D[i] = data
    catch e
        display(e)
    end
end

initdist = findpeak1d(D[22].s.f0,D[22].s.ndisk; eps=D[22].s.eps,tand=D[22].s.tand,granularity=10_000,deviation=0.3)
d0 = initdist*ones(D[22].s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);

b = best(data);
histogram(-data.obj; xlabel="Minimum Boostfactor β²",title="Distribution of Converged States",
    label=false,xlims=(floor(minimum(-data.obj)),ceil(maximum(-data.obj))))
    
showQuality(data)

showDist(data,100; ndiv=10)

showDistribution(data,d0)

m = mean(data.pos,dims=2)[:]


plot(freqsplot/1e9,boost1d(pos2dist(m),freqsplot; eps=data.s.eps,tand=data.s.tand))



for i in 10:100
    fp = genFreqs(i*1e9+25e6,100e6; n=50)
    # display(plot(fp/1e9,boost1d(best(D[i]).dist[:],fp; eps=D[22].s.eps,tand=D[22].s.tand)))

    no = getPeakNo(D[i],fp)
    display(histogram(no; xlabel="Peak Numbers",title="Frequency: $((i*1e9+25e6)/1e9) GHz",yscale=:log))
end



