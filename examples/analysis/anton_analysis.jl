
using Pkg; Pkg.update()

using Dragoon, Plots

include("tools/tools.jl");

tand = 3e-5
eps = 9.35
f0 = 20.2e9

ndisks = collect(5:15)
bws = collect(50:10:100)
# bws = [50]


B = Dict{Tuple{Int,Int},Vector{Float64}}()

for bw in bws
    for ndisk in ndisks
        data = prepareDataAll1d(getPath(),0; f0=f0,df=bw*1e6,ndisk=ndisk,tand=tand,eps=eps)
        
        b = best(data)

        B[(bw,ndisk)] = b.dist
    end
end

for bw in bws
    freqsplot = genFreqs(f0,2*maximum(bws)*1e6; n=200)

    p = plot(; xlabel="Frequency [GHz]",ylabel="Boostfactor β² × 10³",
        title="Optimized Sapphire Disks $bw MHz")#,legend=false)
    vline!(p,[f0-bw*1e6/2,f0+bw*1e6/2]/1e9; c=:black,label="")

    c = palette([:orange,:blue],length(ndisks))

    for (i,ndisk) in Iterators.reverse(enumerate(ndisks))
        data = prepareDataAll1d(getPath(),0; f0=f0,df=bw*1e6,ndisk=ndisk,tand=tand,eps=eps)

        b = best(data)

        plot!(p,freqsplot/1e9,boost1d(b.dist[:],freqsplot; eps=eps,tand=tand)/1e3,
            label="$(ndisk)",c=c[i])
    end

    display(p)

    savefig(p,"lanthal_$(bw).svg")
end



data = prepareDataAll1d(getPath(),-1500; f0=f0,df=50*1e6,ndisk=15,tand=tand,eps=eps)



initdist = findpeak1d(data.s.f0,data.s.ndisk; eps=data.s.eps,tand=data.s.tand,granularity=20_000,deviation=0.5)
showDistribution(data,initdist*ones(data.s.ndisk); dx=1)
showDist(data,100; ndiv=10)

