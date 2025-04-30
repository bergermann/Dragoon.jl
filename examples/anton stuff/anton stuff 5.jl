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


p1 = plot(freqsplot/1e9,abs.(RB1[:,1]))
p2 = plot(freqsplot/1e9,abs2.(RB1[:,2]))

p3 = plot(freqsplot/1e9,abs.(RB2[:,1]))
p4 = plot(freqsplot/1e9,abs2.(RB2[:,2]))






d = 5e-3*ones(9); d[5] = 4e-3
RB1, RB2 = transfer_matrix_2port(freqsplot,d; eps=eps,tand=tand,thickness=1e-3)


p1 = plot(freqsplot/1e9,abs.(RB1[:,1]))
p2 = plot(freqsplot/1e9,abs2.(RB1[:,2]))

p3 = plot(freqsplot/1e9,abs.(RB2[:,1]))
p4 = plot(freqsplot/1e9,abs2.(RB2[:,2]))