
using LinearAlgebra: cross

function tiltinator(z1::Real,z2::Real,z3::Real,r::Real; off::Real=0.)
    θ = deg2rad(off)

    a = [1.5r, sqrt(3)*r/2,z2-z1]
    b = [1.5r,-sqrt(3)*r/2,z3-z1]

    n = cross(b,a)

    tx, ty = [cos(θ) -sin(θ); sin(θ) cos(θ)]*[asin(n[1]/n[3]),asin(n[2]/n[3])]

    return tx, ty
end


a = ones(10)

for i in 1:length(a)
    a[i]
end