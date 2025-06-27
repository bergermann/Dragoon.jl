dx = 0.02
coords = SeedCoordinateSystem(X = -0.5:dx:0.5, Y = -0.5:dx:0.5)

diskR = 0.15

# SetupBoundaries (note that this expects the mirror to be defined explicitly as a region)
epsilon = 24
eps = Array{Complex{Float64}}([NaN, 1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1])
distance = [0.0, 1.00334, 1.0,
                        6.94754, 1.0,
                        7.1766, 1.0,
                        7.22788, 1.0,
                        7.19717, 1.0,
                        7.23776, 1.0,
                        7.07746, 1.0,
                        7.57173, 1.0,
                        7.08019, 1.0,
                        7.24657, 1.0,
                        7.21708, 1.0,
                        7.18317, 1.0,
                        7.13025, 1.0,
                        7.2198, 1.0,
                        7.45585, 1.0,
                        7.39873, 1.0,
                        7.15403, 1.0,
                        7.14252, 1.0,
                        6.83105, 1.0,
                        7.42282, 1.0,
                        0.0]*1e-3

sbdry = SeedSetupBoundaries(coords, diskno=20, distance=distance, epsilon=eps)


# Initialize modes

Mmax = 3
Lmax = 0
modes = SeedModes(coords, ThreeDim=true, Mmax=Mmax, Lmax=Lmax, diskR=diskR)

#  Mode-Vector defining beam shape to be reflected on the system
m_reflect = zeros(Mmax*(2*Lmax+1))
m_reflect[Lmax+1] = 1.0




df = 0.01*1e9
frequencies = 21.98e9:df:22.26e9

# We will build a 3-dim array [reflection / boost factor, mode-vector, frequency ]
# The following function appends to the last dimension
@everywhere zcat(args...) = cat(dims = 3, args...)

# Sweep over frequency
@time for f in frequencies    
    println("Frequency: $f")
    boost, refl = transformer(sbdry,coords,modes; reflect=m_reflect, prop=propagator,diskR=0.15,f=f)
    transpose([boost  refl])
end;
     
