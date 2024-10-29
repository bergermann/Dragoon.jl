





const c0 = 299792458.0


function tm(freqs,distances; tand=0.0,eps=24.0,disc_thickness=1e-3,n_mirror=1e15)
    ϵ = eps*(1.0-1.0im*tand)
    n = sqrt(ϵ)
    A = 1/ϵ
    A0 = 1/sqrt(n_mirror)

    B2 = Vector{Float64}(undef,length(freqs))
    R  = Vector{ComplexF64}(undef,length(freqs))

    
    Gd = [1+n 1-n; 1-n 1+n]/2; Gv = [n+1 n-1; n-1 n+1]/2n
    # G = [1/2*[1+n 1-n; 1-n 1+n], 1/2n*[n+1 n-1; n-1 n+1]]

    S = Matrix{ComplexF64}([(A-1) 0; 0 (A-1)]/2)
    S0 = Matrix{ComplexF64}([(1-A0) 0; 0 (1-A0)]/2)
    M = copy(S)
    G0 = Matrix{ComplexF64}(0.5*[1+n_mirror 1-n_mirror; 1-n_mirror 1+n_mirror])
    T = Matrix{ComplexF64}(0.5*[1+n 1-n; 1-n 1+n])

    W = Matrix{ComplexF64}(undef,2,2)

    for j in eachindex(freqs)
        pd1 = cispi(+2*freqs[j]*disc_thickness)
        pd2 = cispi(-2*freqs[j]*disc_thickness)

        @views for i in Iterators.reverse(eachindex(distances))
            T[:,1] *= pd1
            T[:,2] *= pd2

            mul!(W,T,Gv); T .= W
            # T .= T*Gv

            mul!(W,T,S); M += W
            # M += T*S
        
            T[:,1] *= cispi(+2*freqs[j]*n*distances[i])
            T[:,2] *= cispi(-2*freqs[j]*n*distances[i])

            mul!(W,T,Gd); T .= W
            # T .= T*Gd

            if i > 1
                mul!(W,T,-S); M += W
                # M += T*(-S)
            else
                mul!(W,T,S0); M += W
                # M += T*S0
            end
        end
        
        mul!(W,T,G0); T .= W
        # T .= T*G0

        R[j]  = T[1,2]/T[2,2]
        B2[j] = abs2(M[1,1]+M[1,2])

        # display(T)
        # display(M)

        M .= S
        T .= 1; T[1,1] += n; T[2,2] += n; T[2,1] -= n; T[1,2] -= n; T *= 0.5
    end

    return R, B2
end

r, b = tm(freqs,d; n_mirror=1e30)
