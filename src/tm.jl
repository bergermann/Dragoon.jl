

include("Analytical1D.jl")
using LinearAlgebra: mul!, det
using StaticArrays
using BenchmarkTools
using Plots


const c0 = 299792458.

d = [1.00334, 6.94754, 7.1766, 7.22788, 7.19717,
    7.23776, 7.07746, 7.57173, 7.08019, 7.24657,
    7.21708, 7.18317, 7.13025, 7.2198, 7.45585,
    7.39873, 7.15403, 7.14252, 6.83105, 7.42282]*1e-3
d_ = [d; 0]

freqs = collect(range(21.95,22.1,100))*1e9

ref0 = disk_system(freqs; num_disk=20,tand=0.,disk_epsilon=24.0,disk_thickness=1e-3,spacings=d_)[1]
@btime disk_system($freqs; num_disk=20,tand=0.,disk_epsilon=24.0,disk_thickness=1e-3,spacings=$d_);


function tm(freqs,distances; tand=0.0,eps=24.0,disc_thickness=1e-3,nm=1e15)
    ϵ  = eps*(1.0-1.0im*tand)
    nd = sqrt(ϵ); nm = complex(nm)
    ϵm = nm^2
    A  = 1-1/ϵ
    A0 = 1-1/ϵm

    l1 = length(freqs)
    RB  = MMatrix{2,l1,ComplexF64}(undef)

    Gd = SMatrix{2,2,ComplexF64}((1+nd)/2,   (1-nd)/2,   (1-nd)/2,   (1+nd)/2)
    Gv = SMatrix{2,2,ComplexF64}((nd+1)/2nd, (nd-1)/2nd, (nd-1)/2nd, (nd+1)/2nd)
    G0 = SMatrix{2,2,ComplexF64}((1+nm)/2,   (1-nm)/2,   (1-nm)/2,   (1+nm)/2)
    T  = MMatrix{2,2,ComplexF64}(undef); T .= Gd

    S  = SMatrix{2,2,ComplexF64}( A/2, 0.0im, 0.0im,  A/2)
    S0 = SMatrix{2,2,ComplexF64}(A0/2, 0.0im, 0.0im, A0/2)
    M  = MMatrix{2,2,ComplexF64}(undef); M .= S

    W  = MMatrix{2,2,ComplexF64}(undef) # work matrix for multiplication

    for j in eachindex(freqs)
        pd1 = cispi(-2*freqs[j]*nd*disc_thickness/c0)
        pd2 = cispi(+2*freqs[j]*nd*disc_thickness/c0)

        for i in Iterators.reverse(eachindex(distances))
            T[:,1] .*= pd1
            T[:,2] .*= pd2 # T = Gd*Pd

            mul!(W,T,S); M .-= W    # M = Gd*Pd*S_-1

            mul!(W,T,Gv); T .= W    # T *= Gd*Pd*Gv

            T[:,1] .*= cispi(-2*freqs[j]*distances[i]/c0)
            T[:,2] .*= cispi(+2*freqs[j]*distances[i]/c0)   # T = Gd*Pd*Gv*Gd*S_-1

            if i > 1
                mul!(W,T,S); M .+= W
                mul!(W,T,Gd); T .= W
            else
                mul!(W,T,S0); M .+= W
                mul!(W,T,G0); T .= W
            end
        end
        
        RB[j] = T[1,2]/T[2,2]
        RB[l1+j] = M[1,1]+M[1,2]-(M[2,1]+M[2,2])*T[1,2]/T[2,2]

        # display(det(T))

        M .= S
        T .= 1.0+0.0im; T[1,1] += nd; T[2,2] += nd; T[2,1] -= nd; T[1,2] -= nd; T .*= 0.5
    end

    return RB
end

rb = @btime tm($freqs,$d; nm=1e20);
rb = tm(freqs[1],d; nm=1e20); ref = rb[1,:]; boost = rb[2,:];

plot(freqs/1e9,abs2.(boost))
# plot(freqs/1e9,abs.(ref-ref0))
