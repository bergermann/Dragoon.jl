using Plots
using BoostFractor
using FunctionZeros, SpecialFunctions



coords = SeedCoordinateSystem()
modes = SeedModes(coords; ThreeDim=true,Mmax=3,Lmax=3)

modes.mode_patterns

for m in 1:3, l in 1:7
    showField(modes.mode_patterns[m,l,:,:,1]; title="M=$m, L=$(l-4)")
end

E0 = similar(modes.mode_patterns[1,1,:,:,1]).*0 .+1
E0 .*= [x^2+y^2 <= 0.15^2 for x in coords.X, y in coords.Y]

showField(E0)

coeffs = field2modes(E0,coords,modes)

E1 = modes2field(coeffs,coords,modes)

showField(E1)


Na = sum(abs2.(E0))
mask = [x^2+y^2<0.15^2 for x in coords.X, y in coords.Y]


function mode(m::Integer,l::Integer,coords::Coordinates)
    kr = besselj_zero(l,m)/coords.diskR

    mode = @. besselj(l,kr*coords.R)*cis(-l*coords.Φ)
    @. mode *= coords.diskmaskin
    mode ./= sqrt(sum(abs2.(mode)))
    mode = reshape(mode,(size(mode)...,1))

    return kr, mode
end

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
    coeffs[m,l+n+1] = sum(@. conj(M[m,l+n+1,:])*e0)/na
end

e1 = zeros(ComplexF64,length(x))
for m in 1:n, l in -n:n
    e1 += na*coeffs[m,l+n+1]*M[m,l+n+1,:]
end
vline!(plot(x,abs.(e1);),[-diskR,diskR]); plot!(x,abs.(e0))


plot(x,abs.(Y))
plot(x,abs.(Y); ylim=[-0.05,1.0])
vline!([-diskR,diskR])



a = besselj.(1,x)
b = besselj.(2,x)
