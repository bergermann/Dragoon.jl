
using Pkg; Pkg.update()

using Dragoon
using JLD2, Peaks

include("tools/tools.jl");
@load "examples\\full_20_24.0_6.0e-5.jld2"


D = Dict{Int,Data}()

T = []

for i in 10:27
    try
        threshold = B0[i-9]*0.8
        data = prepareDataAll1d(getPath(),threshold; f0=i*1e9+25e6,df=50e6,tand=6e-5);

        # D[i] = data
        append!(T,data.opttime)
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



f = 100
initdist = findpeak1d(D[f].s.f0,D[f].s.ndisk; eps=D[f].s.eps,tand=D[f].s.tand,granularity=20_000,deviation=0.5)
showDistribution(D[f],initdist*ones(D[f].s.ndisk); dx=0.5)
showDist(D[f],100; ndiv=10)



m = mean(data.pos,dims=2)[:]


plot(freqsplot/1e9,boost1d(pos2dist(m),freqsplot; eps=data.s.eps,tand=data.s.tand))



for i in 10:19
    fp = genFreqs(i*1e9+25e6,100e6; n=50)
    # display(plot(fp/1e9,boost1d(best(D[i]).dist[:],fp; eps=D[22].s.eps,tand=D[22].s.tand)))

    no = getPeakNo(D[i],fp; threshold=0.2)
    # display(histogram(no; xlabel="Peak Numbers",title="Frequency: $((i*1e9+25e6)/1e9) GHz",yscale=:log))
    if any(no.==4); println("quad peak detected"); end
    r = sum(no.==3)
    println("f = $i GHz, ratio = $r")
end

d_ = D[19][findall(no.==3)]
b_ = best(D[19])

fp = genFreqs(19.025e9,150e6; n=1000);
plot(fp/1e9,boost1d(d_.dist[:],fp; eps=data.s.eps,tand=data.s.tand))
plot!(fp/1e9,boost1d(b_.dist[:],fp; eps=data.s.eps,tand=data.s.tand))



wiggle(d_.pos[:],1e-6,1_000,fp,(19e9,19.05e9); eps=data.s.eps,tand=data.s.tand);
wiggle(b_.pos[:],1e-6,1_000,fp,(19e9,19.05e9); eps=data.s.eps,tand=data.s.tand);
wiggle(d_.pos[:],5e-6,1_000,fp,(19e9,19.05e9); eps=data.s.eps,tand=data.s.tand);
wiggle(b_.pos[:],5e-6,1_000,fp,(19e9,19.05e9); eps=data.s.eps,tand=data.s.tand);

wigglewiggle(d_.pos[:],collect(1:10)*1e-6,100,fp,(19e9,19.05e9); eps=data.s.eps,tand=data.s.tand)
wigglewiggle(b_.pos[:],collect(1:10)*1e-6,100,fp,(19e9,19.05e9); eps=data.s.eps,tand=data.s.tand)


