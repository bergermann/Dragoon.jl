
using Pkg; Pkg.update()

using Dragoon, Plots

include("tools/tools.jl");

tand = 6e-5
eps = 24.0
f0 = 22.025e9

ndisks = collect(15:40)
# bws = collect(50:10:100)
bw = 50


B = Dict{Int,Vector{Float64}}()

for ndisk in ndisks
    data = prepareDataAll1d(getPath(),0; f0=f0,df=bw*1e6,ndisk=ndisk,tand=tand,eps=eps)
    
    b = best(data)

    B[ndisk] = b.dist
end





freqsplot = genFreqs(f0,2*maximum(bw)*1e6; n=200)

p = plot(; xlabel="Frequency [GHz]",ylabel="Boostfactor β² × 10³",
    title="Optimized Sapphire Disks $bw MHz")#,ylim=[-1,50])#,legend=false)
vline!(p,[f0-bw*1e6/2,f0+bw*1e6/2]/1e9; c=:black,linestyle=:dash,label="")

c = palette([:orange,:blue],length(ndisks))

for (i,ndisk) in Iterators.reverse(enumerate(ndisks))
    plot!(p,freqsplot/1e9,boost1d(B[ndisk],freqsplot; eps=eps,tand=tand)/1e3,
        label=ndisk%5==0 ? "$(ndisk)" : "",c=c[i])
end

display(p)

savefig(p,"lanthal_$(bw).svg")



wiggle(dist2pos(B[20]),1e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
# wiggle(dist2pos(B[25]),1e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
# wiggle(dist2pos(B[30]),1e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
# wiggle(dist2pos(B[40]),1e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);

wiggle(dist2pos(B[20]),5e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
# wiggle(dist2pos(B[25]),5e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
# wiggle(dist2pos(B[30]),5e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
# wiggle(dist2pos(B[40]),5e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);

wigglewiggle(dist2pos(B[20]),collect(1:10)*1e-6,1000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand)
wigglewiggle(dist2pos(B[25]),collect(1:10)*1e-6,1000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand)
wigglewiggle(dist2pos(B[30]),collect(1:10)*1e-6,1000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand)
wigglewiggle(dist2pos(B[35]),collect(1:10)*1e-6,1000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand)
wigglewiggle(dist2pos(B[40]),collect(1:10)*1e-6,1000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand)



data = prepareDataAll1d(getPath(),-10_000; f0=f0,df=bw*1e6,ndisk=30,tand=tand,eps=eps)



initdist = findpeak1d(data.s.f0,data.s.ndisk; eps=data.s.eps,tand=data.s.tand,granularity=20_000,deviation=0.3)
showDistribution(data,initdist*ones(data.s.ndisk); dx=1)
showDist(data,100; ndiv=10)




data = prepareDataAll1d(getPath(),-10_000; f0=f0,df=bw*1e6,ndisk=40,tand=tand,eps=eps)

fp = genFreqs(f0,2*maximum(bw)*1e6; n=50)
no = Dragoon.getPeakNo(data,fp)
sum(no.==1)
d_ = best(data[findall(no.==3)])




wiggle(d_.pos[:],1e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
wiggle(d_.pos[:],5e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
wiggle(b_.pos[:],1e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);
wiggle(b_.pos[:],5e-6,10_000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand);

wigglewiggle(d_.pos[:],collect(1:10)*1e-6,1000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand)
wigglewiggle(b_.pos[:],collect(1:10)*1e-6,1000,freqsplot,(22e9,22.05e9); eps=eps,tand=tand)

