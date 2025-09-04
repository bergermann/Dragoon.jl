
using BoostFractor



dx = 0.01
coords = SeedCoordinateSystem(X = -0.5:dx:0.5, Y = -0.5:dx:0.5)

diskR = 0.15

# SetupBoundaries
eps = Array{Complex{Float64}}([NaN,1])
distance = [0.0, 1.5, ]*1e-3
tilt_x = [0.,0*deg2rad(0.1)]
tilt_y = [0.,0.]

sbdry = SeedSetupBoundaries(coords, diskno=0, distance=distance, epsilon=eps,
    relative_tilt_x=tilt_x,relative_tilt_y=tilt_y)

# Initialize modes
Mmax = 4
Lmax = 3 # l-Modes are irrelevant for the azimuthally symmetric haloscope
# For a 1D calculation:
#modes = SeedModes(coords, ThreeDim=false, Mmax=Mmax, Lmax=Lmax, diskR=diskR)
# For 3D:
modes = SeedModes(coords, ThreeDim=true, Mmax=Mmax, Lmax=Lmax, diskR=diskR)


freqs = collect(range(1,30,1000)*1e9)

b = zeros(ComplexF64,Mmax*(2Lmax+1),length(freqs))

@time for i in eachindex(freqs)
    b[:,i] = transformer(sbdry,coords,modes; prop=propagator,f=freqs[i],diskR=diskR)
end

B = sum(abs2,b; dims=1)[:]

plot(freqs,B)