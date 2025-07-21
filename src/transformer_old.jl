

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

# s = setup(d,4,3; dx=0.02)

# tilt!(s.sbdry,0.001)

# B, R, b, b_sum, r = boost(freqs,s)


# display(plot(freqs/1e9,b; title="transformer",legend=false))
# display(plot(freqs/1e9,b_sum; title="transformer",legend=false))

function setup_prop(freqs,distances,tilts,(coords,sbdry,modes,m_reflect,diskR))
    pg = calc_propagation_matrices_grid(sbdry,coords,modes,distances,freqs;
        diskR=diskR,tilt_x_grid=deg2rad.(tilts),tilt_y_grid=deg2rad.(tilts));

    return [pg[i,j,1,1,1,:,:] for i in eachindex(sbdry.distance), j in eachindex(freqs)]
end

function boost(freqs,distances,tiltsx,tiltsy,(coords,sbdry,modes,m_reflect,diskR),p)
    tilt!(sbdry,tiltsx,tiltsy)
    move!(sbdry,distances)

    
    prop_matrices_set_interp = interpolate_prop_matrix(itp,dist_shift);

    Eout = calc_boostfactor_modes(sbdry,coords,modes,freqs,p;diskR=diskR);
    
    return Eout 
end



s = setup(d,2,1; dx=0.02)
@time p = setup_prop(freqs,0,-0.005:0.001:0.005,s);

boost(freqs,d,zeros(20),zeros(20),s,p)