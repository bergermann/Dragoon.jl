

using BoostFractor, Plots, BenchmarkTools

include("transformer_optim_utilities.jl")
include("transformer_utils.jl")



freqs = collect(range(21.98e9,22.26e9,30)); λ = 299792458.0/22e9
d = [
    1.00334,
    6.94754,
    7.17660,
    7.22788,
    7.19717,
    7.23776,
    7.07746,
    7.57173,
    7.08019,
    7.24657,
    7.21708,
    7.18317,
    7.13025,
    7.21980,
    7.45585,
    7.39873,
    7.15403,
    7.14252,
    6.83105,
    7.42282
]*1e-3;

s = setup(d,4,3; dx=0.02)

tilt!(s.sbdry,0.001)

B, R, b, b_sum, r = boost(freqs,s)


display(plot(freqs/1e9,b; title="transformer",legend=false))
display(plot(freqs/1e9,b_sum; title="transformer",legend=false))




calc_propagation_matrices_grid(sbdry,coords,modes,spacing_grid,frequencies;tilt_x_grid=0,tilt_y_grid=0, prop=propagator, diskR=0.15)


prop_matrix_grid = calc_propagation_matrices_grid(sbdry,coords,modes,0,frequencies;prop=prop, diskR=diskR);
prop_matrix = [prop_matrix_grid[r,f,1,1,1,:,:] for r=1:n_region, f=1:length(frequencies)]
Eout_init = calc_boostfactor_modes(sbdry,coords,modes,frequencies, prop_matrix;prop=prop, diskR=diskR);


# @btime Eout_init = calc_boostfactor_modes($sbdry,$coords,$modes,$frequencies,$prop_matrix; prop=prop, diskR=diskR);


# Eout_init