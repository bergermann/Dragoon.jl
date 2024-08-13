

function scanSingleDiscs(d0::Vector{Float64},f0::Float64,df::Float64,nfreqs::Int,
        n1::Int,n2::Int)

    freqs = genFreqs(f0,df; n=nfreqs)
    l = λ(f0)

    for j in eachindex(d0)
        O = zeros(n1)
        d_ = copy(d0); d_[j] = 0
        d = zeros(length(d0)); d[j] = 1

        for i in eachindex(O)
            O[i] = minimum(boost1d(d_+d*l*i/n2,freqs))
        end

        display(plot((1:n1)/n2,O; xlabel="d_$j/λ",ylabel="minimum boost"))
    end

    return
end

