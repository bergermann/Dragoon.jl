


using LinearAlgebra: mul!, det
using StaticArrays


const c0 = 299792458.


function tm(freqs,distances; tand=0.0,eps=24.0,disc_thickness=1e-3,nm=1e15)
    ϵ  = eps*(1.0-1.0im*tand)
    nd = sqrt(ϵ)
    A  = 1/ϵ-1
    A0 = complex(1-1/sqrt(nm))

    # RB  = Matrix{ComplexF64}(undef,2,length(freqs))

    # Gd = ComplexF64[1+nd 1-nd; 1-nd 1+nd]; Gd ./= 2
    # Gv = ComplexF64[nd+1 nd-1; nd-1 nd+1]; Gv ./= 2nd
    # G0 = ComplexF64[1+nm 1-nm; 1-nm 1+nm]; G0 ./= 2
    # T  = ComplexF64[1+nd 1-nd; 1-nd 1+nd]; T  ./= 2

    # S  = ComplexF64[ A 0.0im; 0.0im  A]; S  ./= 2
    # S0 = ComplexF64[A0 0.0im; 0.0im A0]; S0 ./= 2
    # M  = copy(S)

    # W  = Matrix{ComplexF64}(undef,2,2)

    l1 = length(freqs); l2 = 2*l1
    RB  = MMatrix{2,l1,ComplexF64}(undef)

    Gd = SMatrix{2,2,ComplexF64}((1+nd)/2, (1-nd)/2, (1-nd)/2, (1+nd)/2) # ; Gd ./= 2
    Gv = SMatrix{2,2,ComplexF64}((nd+1)/2nd, (nd-1)/2nd, (nd-1)/2nd, (nd+1)/2nd) # ; Gv ./= 2nd
    G0 = SMatrix{2,2,ComplexF64}((1+nm)/2, (1-nm)/2, (1-nm)/2, (1+nm)/2) # ; G0 ./= 2
    T  = MMatrix{2,2,ComplexF64}(undef); T .= Gd

    S  = SMatrix{2,2,ComplexF64}( A/2, 0.0im, 0.0im,  A/2) # ; S  ./= 2
    S0 = SMatrix{2,2,ComplexF64}(A0/2, 0.0im, 0.0im, A0/2) # ; S0 ./= 2
    M  = MMatrix{2,2,ComplexF64}(undef); M .= S

    W  = MMatrix{2,2,ComplexF64}(undef) # work matrix for multiplication

    d = 0

    @views for j in eachindex(freqs)
        pd1 = cispi(+2*freqs[j]*nd*disc_thickness/c0)
        pd2 = cispi(-2*freqs[j]*nd*disc_thickness/c0)

        for i in Iterators.reverse(eachindex(distances))
            T[:,1] .*= pd1
            T[:,2] .*= pd2

            d += disc_thickness

            display(d)

            mul!(W,T,Gv); T .= W

            mul!(W,T,S); M .+= W
        
            T[:,1] .*= cispi(+2*freqs[j]*distances[i]/c0)
            T[:,2] .*= cispi(-2*freqs[j]*distances[i]/c0)

            d += distances[i]
            
            display(d)

            if i > 1
                mul!(W,T,Gd); T .= W
                mul!(W,T,S); M .-= W
            else
                mul!(W,T,G0); T .= W
                mul!(W,T,S0); M .+= W
            end
        end
        
        RB[1,j] = T[1,2]/T[2,2]
        RB[2,j] = M[1,1]+M[1,2]

        # display(det(T))

        M .= S
        T .= 1.0+0.0im; T[1,1] += nd; T[2,2] += nd; T[2,1] -= nd; T[1,2] -= nd; T .*= 0.5
    end

    return RB
end

rb = @btime tm($freqs,$d; nm=1e10)
rb = tm(freqs[1],d; nm=1e9); ref = rb[1,:]; boost = rb[2,:];

plot(freqs/1e9,abs2.(boost))
plot(freqs/1e9,abs.(ref-ref0))


function test(f)
    l1 = length(f)

    MMatrix{2,l1,ComplexF64}(undef)

    return
end

@btime test($freqs)
test(freqs)