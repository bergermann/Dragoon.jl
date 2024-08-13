

function wiggle(pos::Vector{Float64}, sigx::Real, n::Int,
            freqsplot::Vector{Float64}, bounds::Tuple{Float64,Float64};
            showplots::Bool=true, ndiv::Int=10)

        b0 = boost1d(pos2dist(pos), freqsplot)
        obj0 = minimum(b0[@.(bounds[1] < freqsplot < bounds[2])])

        db = zeros(Float64, length(freqsplot), n)
        obj = zeros(Float64, n)
        dsum = zeros(Float64, n)
        dsumabs = zeros(Float64, n)
        dmax = zeros(Float64, n)

        R = Vector{Float64}[]

        for i in 1:n
            r = sigx * randn(length(pos))
            dsum[i] = sum(r)
            dsumabs[i] = sum(abs.(r))
            dmax[i] = maximum(abs.(r))

            db[:, i] = boost1d(pos2dist(pos + r), freqsplot)
            obj[i] = minimum(db[@.(bounds[1] < freqsplot < bounds[2]), i])

            if obj[i] < obj0
                push!(R, pos2dist(pos + r))
            end
        end

        db .-= b0

        q = collect(0:1/ndiv:1)
        Q = zeros(length(q), length(freqsplot))
        for i in axes(Q, 2)
            Q[:, i] = quantile(db[i, :], q)
        end

        if showplots
            p1 = plot(; xlabel="Frequency [GHz]", ylabel="Boostfactor β²")

            for i in 1:div(length(q) - 1, 2)
                plot!(p1, freqsplot / 1e9, b0 + Q[i, :]; linewidth=0,
                    fillrange=b0 + Q[end-i+1, :], fillalpha=0.3, fillcolor=1, label=""
                )
            end

            plot!(p1, freqsplot / 1e9, b0; label="original")
            plot!(p1, freqsplot / 1e9, b0 + Q[div(length(q) - 1, 2), :], label="50% quant")

            p2 = histogram(obj)
            vline!(p2, [obj0])

            p3 = histogram2d(obj, dsum / 1e-6; bins=20, xlabel="objective value", ylabel="∑Δp_i [mm]")
            p4 = histogram2d(obj, dsumabs / 1e-6; bins=20, xlabel="objective value", ylabel="∑|Δp_i| [mm]")
            p5 = histogram2d(obj, dmax / 1e-6; bins=20, xlabel="objective value", ylabel="max(Δp_i) [mm]")

            display(p1)
            display(p2)
            display(p3)
            display(p4)
            display(p5)
        end

        return Q, R, obj
end

function wigglewiggle(pos::Vector{Float64},sigxs::Vector{<:Real},n::Int,
    freqsplot::Vector{Float64},bounds::Tuple{Float64,Float64}; ndiv::Int=10)

    q = collect(0:1/ndiv:1)
    Q = zeros(Float64,length(sigxs),length(q))

    for i in eachindex(sigxs)
        println("Step $i/$(length(sigxs))")
        _, _, obj = wiggle(pos,sigxs[i],n,freqsplot,bounds; showplots=false,ndiv=ndiv)

        Q[i,:] = quantile(obj,q)
    end

    p1 = plot(; xlabel="σ_Δp [μm]",ylabel="Boostfactor β²",legend=false,
        title="Boostfactor Deteration on Positional Uncertainty")

    for i in axes(Q,2)
        plot!(p1,sigxs/1e-6,Q[:,i]; c=1,
            lw=0.5,fillrange=Q[:,end],fillcolor=1,fillalpha=0.2)
    end

    display(p1)

    return
end
