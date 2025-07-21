

using BoostFractor

function setup(distances,M,L; eps=24,Δx=1.,dx=0.02,diskR=0.15)
    ndisk = length(distances)

    X = collect(-Δx/2:dx:Δx/2)
    coords = SeedCoordinateSystem(X=X,Y=copy(X))

    d = zeros(Float64,2ndisk+2); d[2:2:end-2] .= distances; d[3:2:end-1] .= 1e-3
    e = zeros(ComplexF64,2ndisk+2); e[2:2:end] .= 1.; e[3:2:end] .= eps; e[1] = NaN

    sbdry = SeedSetupBoundaries(coords,diskno=ndisk,distance=d,epsilon=e)
    modes = SeedModes(coords,ThreeDim=true,Mmax=M,Lmax=L,diskR=diskR)

    m_reflect = zeros(Float64,M*(2*L+1))
    m_reflect[L+1] = 1.0

    return (coords=coords, sbdry=sbdry, modes=modes, m_reflect=m_reflect, disk=diskR)
end

function tilt!(sbdry::SetupBoundaries,deg::Real)
    maxtilt = deg2rad(deg)

    for i in 1:20
        sbdry.relative_tilt_x[2*i] = sbdry.relative_tilt_x[1+2*i] = maxtilt*(2*rand()-1)
        sbdry.relative_tilt_y[2*i] = sbdry.relative_tilt_y[1+2*i] = maxtilt*(2*rand()-1)
    end

    return
end

function tilt!(sbdry::SetupBoundaries,tiltsx::Vector{<:Real},tiltsy::Vector{<:Real})
    @assert length(tiltsx) == length(tiltsy) == Int((length(sbdry.distance)-2)/2)
        "Tilt arrays needs entry for each disc."

    for i in 1:20
        sbdry.relative_tilt_x[2*i] = sbdry.relative_tilt_x[1+2*i] = deg2rad(tiltsx[i])
        sbdry.relative_tilt_y[2*i] = sbdry.relative_tilt_y[1+2*i] = deg2rad(tiltsy[i])
    end

    return
end

function move!(sbdry::SetupBoundaries,distances::Vector{<:Real})
    @assert length(distances) == Int((length(sbdry.distance)-2)/2)
        "Distance array needs entry for each disc."

    sbdry.distance[2:2:end-2] .= distances

    return
end

function boost(freqs,(coords,sbdry,modes,m_reflect,diskR))
    ML = modes.M*(2modes.L+1)
    B = Vector{Vector{ComplexF64}}(undef,length(freqs))
    R = Vector{Vector{ComplexF64}}(undef,length(freqs))

    @time for i in eachindex(freqs)
        B[i], R[i] = transformer(sbdry,coords,modes; reflect=m_reflect,diskR=diskR,f=freqs[i])
    end

    b = [[abs2(B[i][j]) for i in eachindex(B)] for j in 1:ML]
    b_sum = [sum(abs2,b) for b in B]
    
    r = [[R[i][j] for i in eachindex(R)] for j in 1:ML]

    return B, R, b, b_sum, r
end



# freqs = collect(range(21.98e9,22.26e9,10)); λ = 299792458.0/22e9
# d = [
#     1.00334,
#     6.94754,
#     7.17660,
#     7.22788,
#     7.19717,
#     7.23776,
#     7.07746,
#     7.57173,
#     7.08019,
#     7.24657,
#     7.21708,
#     7.18317,
#     7.13025,
#     7.21980,
#     7.45585,
#     7.39873,
#     7.15403,
#     7.14252,
#     6.83105,
#     7.42282
# ]*1e-3;

# s = setup(d,4,3; dx=0.02)

# tilt!(s.sbdry,0.005)

# B, R, b, b_sum, r = boost(freqs,s)


# plot(freqs/1e9,b; title="transformer",legend=false)
# plot(freqs/1e9,b_sum; title="transformer",legend=false)

