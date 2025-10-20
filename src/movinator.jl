

function movinator(old::Vector{<:Real},new::Vector{<:Real};
        thickness::Real=0.1,lbl::Real=-Inf64,ubr::Real=Inf64)

    l = length(old)

    @assert l == length(new) "Old and new positions need same lengths."
    @assert issorted(old) && issorted(new) "Both position vectors need to be in ascending order."

    old_ = copy(old); 
    idx = zeros(Int,l)

    for j in eachindex(idx)
        for i in eachindex(idx)
            if idx[i] != 0; continue; end

            lb = i == 1 ? lbl : old_[i-1]+thickness
            ub = i == l ? ubr : old_[i+1]

            if lb < new[i] < ub
                old_[i] = new[i]
                idx[i] = j

                break
            end
        end
    end

    @assert all(idx .!= 0) "Operation unsuccessfull. Not everything moved."
    @assert allunique(idx) "Operation unsuccessfull. No unique steps."

    return idx
end

function genpos(n::Int,dmax::Real=0.1; thickness::Real=0.1)
    @assert n > 0 "Disc number needs to be positive integer."
    @assert dmax > 0 "dmax needs to be positive."
    @assert thickness > 0 "thickness needs to be positive."

    p = zeros(Float64,n)

    p[1] = rand()*dmax

    for i in 2:n
        p[i] = p[i-1]+thickness+rand()*dmax
    end

    return p
end

a = genpos(10)
b = genpos(10,0.2)

movinator(a,b)

for i in 1:100_000
    a = genpos(100)
    b = genpos(100,0.2)
    
    movinator(a,b)
end