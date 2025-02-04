

using Plots
using BoostFractor
using FunctionZeros, SpecialFunctions




n = 10
coords = SeedCoordinateSystem(; X=-0.5:0.01:0.5,Y=-0.5:0.01:0.5)
modes = SeedModes(coords; ThreeDim=true,Mmax=n,Lmax=n)
mask = [x^2+y^2<0.15^2 for x in coords.X, y in coords.Y];

# modes.mode_patterns

# for m in 1:n l in 1:2n+1
#     modes.mode_patterns[m,l,:,:,1] .*= mask
# end
    


for m in 1:3, l in 1:2n+1
    showCross(modes.mode_patterns[m,l,:,:,1]; title="M=$m, L=$(l-4)")
end



E0 = similar(modes.mode_patterns[1,1,:,:,1]).*0 .+1
E0 .*= mask
Na = sqrt(sum(abs2.(E0))); #E0 ./= Na

# showCross(E0)



coeffs1 = field2modes(E0,coords,modes); sum(abs2.(coeffs1))
E1 = modes2field(coeffs1,coords,modes); showCross(E1); plot!(abs.(E0[501,:]))

coeffs = zeros(ComplexF64,n,2n+1)
for m in 1:n, l in 1:2n+1
    coeffs[m,l] = sum(conj.(modes.mode_patterns[m,l,:,:,1]).*E0)
end; coeffs



E2 = similar(modes.mode_patterns[1,1,:,:,1]).*0
for m in 1:n, l in 1:2n+1
    E2 += coeffs[m,l]*modes.mode_patterns[m,l,:,:,1]
end; showCross(E2)

# showField(E1)











function showCross(E::Array{ComplexF64}; kwargs...)
    n = size(E,1); @assert isodd(n) "Grid for E needs odd edge lengths."
    n2 = div(n+1,2)

    for i in axes(E,3)
        display(plot(abs.(E[n2,:,i]); kwargs...))
    end
    
    return
end








# function mode(m::Integer,l::Integer,coords::Coordinates)
#     kr = besselj_zero(l,m)/coords.diskR

#     mode = @. besselj(l,kr*coords.R)*cis(-l*coords.Φ)
#     @. mode *= coords.diskmaskin
#     mode ./= sqrt(sum(abs2.(mode)))
#     mode = reshape(mode,(size(mode)...,1))

#     return kr, mode
# end

p1 = modes.mode_patterns[1,1,:,:,1].*mask
p2 = modes.mode_patterns[3,5,:,:,1].*mask

showField(p1)

sum(abs2.(p1))
sum(abs2.(p2))
sum(@. conj(p1)*p2)


n = 10
x = collect(-0.5:0.001:0.5); diskR = 0.15
Y = zeros(ComplexF64,length(x))
Φ = zeros(Float64,length(x)); Φ[x.<0] .= -π
mask = @. abs(x) <= diskR
e0 = ones(ComplexF64,size(x)).*mask; na = sum(abs2.(e0))
M = zeros(ComplexF64,n,2n+1,length(x))
for m in 1:n, l in -n:n
    kr = besselj_zero(l,m)/diskR
    y = @. besselj(l,kr*x)*cis(-l*Φ)*mask
    y /= sum(abs2.(y))
    M[m,l+n+1,:] .= y
    # display(vline!(plot(x,abs.(y); ylim=[-0.05,1.0]),[-diskR,diskR]))
end

coeffs = zeros(ComplexF64,n,2n+1)
for m in 1:n, l in -n:n
    coeffs[m,l+n+1] = (sum(@. conj(M[m,l+n+1,:])*e0)/na)^2
end

e1 = zeros(ComplexF64,length(x))
for m in 1:n, l in -n:n
    e1 += na^2*coeffs[m,l+n+1]*M[m,l+n+1,:]
end
vline!(plot(x,abs.(e1);),[-diskR,diskR]); plot!(x,abs.(e0))


plot(x,abs.(Y))
plot(x,abs.(Y); ylim=[-0.05,1.0])
vline!([-diskR,diskR])



a = besselj.(1,x)
b = besselj.(2,x)
