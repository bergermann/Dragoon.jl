


a = rand(10)
b = rand(10)*0.5 .+ 0.5

sort!(a); sort!(b)

function movinator(old::Vector{<:Real},new::Vector{<:Real})
    l = length(old)

    @assert l == length(new) "Old and new positions need same lengths."
    @assert issorted(old) && issorted(new) "Both position vectors need to be in ascending order."

    old_ = copy(old); 
    idx = zeros(Int,l)

    for j in eachindex(idx)
        for i in eachindex(idx)
            if idx[i] != 0; continue; end

            lb = i == 1 ? -Inf64 : old_[i-1]
            ub = i == l ?  Inf64 : old_[i+1]

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

for i in 1:100_000
    a = rand(10)
    x = rand()
    b = rand(10)*x .+ (1-x)

    sort!(a); sort!(b)
    
    movinator(a,b)
end