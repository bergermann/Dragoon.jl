
using Pkg; Pkg.update()

using Dragoon, Plots

include("tools/tools.jl");

tand = 6e-5
eps = 24.0
f0 = 22.025e9

ndisks = collect(15:30)
# bws = collect(50:10:100)
bws = [50]


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
