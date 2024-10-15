
using Dragoon, Plots

include("tools/tools.jl");

tand = 3e-5
eps = 9.35
f0 = 20.2e9

ndisks = collect(5:15)
bws = collect(50:10:100)*1e6


for bw in bws
    freqsplot = genFreqs(f0,3*bw; n=200)

    p = plot(; xlabel="Frequency [GHz]",ylabel="Boostfactor β² × 10³",
        title="Optimized States 5-15 Sapphire Disks",legend=false)

    for ndisk in ndisks
        data = prepareDataAll1d(getPath(),0; f0=f0,df=bw,tand=tand,eps=eps)

        b = best(data)

        plot!(p,freqsplot/1e9,boost1d(b.dist[:],freqsplot; eps=eps,tand=tand)/1e3,
            label="")
    end

    break
end
