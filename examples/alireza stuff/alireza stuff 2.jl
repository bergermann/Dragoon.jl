using Pkg

# Pkg.add(url="https://github.com/mppmu/BoostFractor.jl.git")
# Pkg.add(url="https://github.com/bergermann/Dragoon.jl.git")
Pkg.update()

using Dragoon
using Plots
using HDF5

eps = 25.
tand = 0
freqsplot = genFreqs([27,33]*1e9; n=1000);



n = 6
d = 5e-3*ones(n);
RB = transfer_matrix(freqsplot,d; eps=eps,tand=tand,thickness=1e-3)

p1 = plot(freqsplot/1e9,abs2.(RB[:,2]); xlabel="Frequency [GHz]",ylabel="Boostfactor β²",legend=false)


h5write("6discs.h5","data/boost",RB[:,2])
h5write("6discs.h5","data/freqs",freqsplot)



n = 10
d = 5e-3*ones(n);
RB = transfer_matrix(freqsplot,d; eps=eps,tand=tand,thickness=1e-3)

p1 = plot(freqsplot/1e9,abs2.(RB[:,2]); xlabel="Frequency [GHz]",ylabel="Boostfactor β²",legend=false)

h5write("10discs.h5","data/boost",RB[:,2])
h5write("10discs.h5","data/freqs",freqsplot)