using Pkg

# Pkg.add(url="https://github.com/mppmu/BoostFractor.jl.git")
# Pkg.add(url="https://github.com/bergermann/Dragoon.jl.git")
Pkg.update()

using Dragoon
using Plots
using HDF5
# using Dates

eps = 9.35
tand = 1e-4
freqsplot = genFreqs([30,36]*1e9; n=1000);


d = 5e-3*ones(9);
RB1, RB2 = transfer_matrix_2port(freqsplot,d; eps=eps,tand=tand,thickness=1e-3)
all(RB1 .== RB2)

p1 = plot(freqsplot/1e9,abs2.(RB1[:,2]); xlabel="Frequency [GHz]",ylabel="Boostfactor β²",legend=false)
p2 = plot(freqsplot/1e9,abs.(RB1[:,1]); xlabel="Frequency [GHz]",ylabel="Reflectivity R",label="abs",c=:red,lw=2)
plot!(p2,freqsplot/1e9,real.(RB1[:,1]); label="real",c=:blue)
plot!(p2,freqsplot/1e9,imag.(RB1[:,1]); label="imag",c=:darkorange)
a = angle.(RB1[:,1])/π .+1;
p3 = plot(freqsplot/1e9,a; xlabel="Frequency [GHz]",ylabel="Phase ∠(R)/π",legend=false)
display(p1); display(p2); display(p3);

savefig(p1,"boost.svg"); savefig(p2,"reflectivity.svg"); savefig(p3,"phase.svg")
h5write("data.h5","data/boost",RB1[:,2])
h5write("data.h5","data/reflectivity",RB1[:,1])
h5write("data.h5","data/freqs",freqsplot)


p4 = plot(freqsplot/1e9,abs2.(RB2[:,2]); xlabel="Frequency [GHz]",ylabel="Boostfactor β²",legend=false)
p5 = plot(freqsplot/1e9,abs.(RB2[:,1]); xlabel="Frequency [GHz]",ylabel="Reflectivity R",label="abs",c=:red,lw=2)
plot!(p5,freqsplot/1e9,real.(RB2[:,1]); label="real",c=:blue)
plot!(p5,freqsplot/1e9,imag.(RB2[:,1]); label="imag",c=:darkorange)
a = angle.(RB2[:,1])/π .+1;
p6 = plot(freqsplot/1e9,a; xlabel="Frequency [GHz]",ylabel="Phase ∠(R)/π",legend=false)
display(p4); display(p5); display(p6);








d = 5e-3*ones(9); d[5] = 4e-3
RB1, RB2 = transfer_matrix_2port(freqsplot,d; eps=eps,tand=tand,thickness=1e-3)
all(RB1 .== RB2)

p1 = plot(freqsplot/1e9,abs2.(RB1[:,2]); xlabel="Frequency [GHz]",ylabel="Boostfactor β²",legend=false)
p2 = plot(freqsplot/1e9,abs.(RB1[:,1]); xlabel="Frequency [GHz]",ylabel="Reflectivity R",label="abs",c=:red,lw=2)
plot!(p2,freqsplot/1e9,real.(RB1[:,1]); label="real",c=:blue)
plot!(p2,freqsplot/1e9,imag.(RB1[:,1]); label="imag",c=:darkorange)
a = angle.(RB1[:,1])/π .+1;
p3 = plot(freqsplot/1e9,a; xlabel="Frequency [GHz]",ylabel="Phase ∠(R)/π",legend=false)
display(p1); display(p2); display(p3);

savefig(p1,"boost.svg"); savefig(p2,"reflectivity.svg"); savefig(p3,"phase.svg")
h5write("data.h5","data/boost",RB1[:,2])
h5write("data.h5","data/reflectivity",RB1[:,1])
h5write("data.h5","data/freqs",freqsplot)
