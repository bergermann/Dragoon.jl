



function showFields(pos::Vector{Float64}, frequency::Float64, freqsplot=nothing;
    R::Float64=0.1, eps::Float64=24.0, tand::Float64=0.0, thickness::Float64=1e-3,
    obj::Float64=Inf64)

    ndisk = 20#length(pos)

    eps_ = ComplexF64[!iseven(i) ? 1 : eps for i in 1:2*ndisk+1]
    prepend!(eps_, NaN)

    d = pos2dist(pos; disk_thickness=thickness)
    distance = Float64[0]
    for i in 1:ndisk
        append!(distance, d[i], thickness)
    end
    append!(distance, 0)

    coords = SeedCoordinateSystem(X=[1e-9], Y=[1e-9])

    sbdry = SeedSetupBoundaries(coords; diskno=ndisk, distance=distance, epsilon=eps_)

    modes = SeedModes(coords, ThreeDim=false, Mmax=1, Lmax=0, diskR=R)
    m_reflect = Float64[1.0]

    boost, refl = transformer(sbdry, coords, modes;
        reflect=m_reflect, prop=propagator1D, diskR=R, f=frequency)
    Eout = transpose([boost refl])

    if isnothing(freqsplot)
        p = plot(; layout=(2, 1), size=(800, 600), ylabel="E/E0", titlelocation=:left,)
    else
        p = plot(; layout=(3, 1), size=(800, 900), ylabel="E/E0", titlelocation=:left,)
    end

    annotate!(p[1], (0.75, -0.125), text(@sprintf("%.3f GHz", frequency / 1e9), :left))
    title!(p[1], "Reflection")
    title!(p[2], "Axion Induced")
    plot!(p[1]; legend=false, xformatter=_ -> "", xlabel="")
    plot!(p[2]; legend=:bottomright, xlabel="z Position [m]")

    drawBdry!(p[1], sbdry)
    drawBdry!(p[2], sbdry)

    full_fields_r = BoostFractor.transformer_trace_back(Eout[2, :], m_reflect, sbdry, coords, modes;
        prop=propagator1D, f=frequency)
    # plot_1d_field_pattern!(p[1],-autorotate(full_fields_r[:,:,1]), sbdry, frequency)
    sr = plot_1d_field_pattern!(p[1], full_fields_r[:, :, 1], sbdry, frequency)


    full_fields_a = BoostFractor.transformer_trace_back(Eout[1, :], Float64[0], sbdry, coords, modes;
        prop=propagator1D, f=frequency, inlcudes_axion=true)
    # plot_1d_field_pattern(-(full_fields_a[:,:,1]), sbdry, frequency,
    #     add_ea=true, overallphase=exp(-1im*pi/2*0.95))
    sa = plot_1d_field_pattern!(p[2], full_fields_a[:, :, 1], sbdry, frequency,
        add_ea=true, overallphase=exp(-1im * pi / 2 * 0.95))

    if !isnothing(freqsplot)
        plot!(p[3], freqsplot / 1e9,
            boost1d(d, freqsplot; eps=eps, tand=tand, thickness=thickness) / 1e3,
            xlabel="Frequency [GHz]", ylabel="Boost β² ×10³", label="", legend=false)
        vline!(p[3], [frequency / 1e9], label="")

        if !isequal(obj, Inf64)
            annotate!(p[3], (0.8, 1.0),
                text(@sprintf("Obj. value:\n%.1f", obj), :left))
        end
    end

    return p, (sr, sa)
end


function plot_1d_field_pattern!(p, full_solution_regions, bdry::SetupBoundaries, f;
    add_ea=false, overallphase=1)

    ztot = 0
    S = 0
    Nregions = length(bdry.eps)
    c0 = 299792458.0

    maxE = maximum(abs.(full_solution_regions[:, :]))
    ylims!(p, Tuple([-1, 1] * max(ylims(p)[2], 2.2 * maxE)))

    for i in 1:Nregions
        kreg = 2pi / c0 * f * sqrt(bdry.eps[i])

        S += plotRegion1d!(p, full_solution_regions[i, 1] .* overallphase,
            full_solution_regions[i, 2] .* overallphase,
            ztot, ztot + bdry.distance[i], kreg;
            extraspace=(i == Nregions),
            Ea=(add_ea ? (1 / bdry.eps[i]) .* overallphase : 0))

        ztot += bdry.distance[i]
    end

    return S
end

function plotRegion1d!(p, R_, L_, z0, z1, k; extraspace=false, Ea=0)
    z1 += 1e-9
    maximum = (z1 + (extraspace ? 10e-3 : 0))
    z = vcat(z0:2e-4:maximum, maximum)
    dz2 = .-(z .- z1)

    Rf = L_ * exp.(+1im .* k .* dz2)
    Lf = R_ * exp.(-1im .* k .* dz2)

    E = Rf .+ Lf .+ Ea

    # plot!(p,z, real.(E); c=c,label=(extraspace ? label : ""))
    plot!(p, z, abs.(E), c=:black, linestyle=:dash, linewidth=0.5,
        label=(extraspace ? "|E|" : ""))
    plot!(p, z, real.(E), c=:blue, label=(extraspace ? "Re(E)" : ""))
    plot!(p, z, imag.(E), c=:red, label=(extraspace ? "Im(E)" : ""))
    # plot!(p,z, angle.(E), c=:black, label=(extraspace ? "∠(E)" : ""))

    return sum(abs.(E))
end

function drawBdry!(p, bdry)
    z0 = 0

    for i in eachindex(bdry.eps)
        z1 = z0 + bdry.distance[i]
        if abs.(bdry.eps[i]) > 100 || bdry.eps[i] == NaN
            vspan!(p, [z0 == 0 ? -0.0025 : z0, z1], c=:darkorange, label="")
        elseif abs.(bdry.eps[i]) != 1
            vspan!(p, [z0 == 0 ? -0.0025 : z0, z1], c=:lightgray, label="")
        end

        z0 = z1
    end

    return
end

function autorotate(full_solution_regions)
    ang = angle.(full_solution_regions[2, 1] .+ full_solution_regions[2, 2])
    #sgn = real.((full_solution_regions[2,1].+full_solution_regions[2,2]).*exp(-1im*ang)) .> 0
    return full_solution_regions .* exp(-1im * ang)
end

function showFields(pos::Vector{Float64}, frequencies::Vector{Float64}, freqsplot=nothing;
    R::Float64=0.1, eps::Float64=24.0, tand::Float64=0.0, thickness::Float64=1e-3,
    yscaling::Bool=true, obj::Float64=Inf64)

    P = [showFields(pos, f, freqsplot; R=R, eps=eps, tand=tand, thickness=thickness, obj=obj)
         for f in frequencies]

    xlim1 = maximum([xlims(p[1][1])[2] for p in P])

    for i in eachindex(P)
        xlims!(P[i][1][1], (-0.002, xlim1))
        xlims!(P[i][1][2], (-0.002, xlim1))
    end

    if yscaling
        ylim1 = maximum([ylims(p[1][1])[2] for p in P])
        ylim2 = maximum([ylims(p[1][2])[2] for p in P])

        for i in eachindex(P)
            ylims!(P[i][1][1], (-ylim1, ylim1))
            ylims!(P[i][1][2], (-ylim2, ylim2))
        end
    end

    for p in P
        display(p[1])
    end

    p_ = plot(; layout=(2, 1), size=(800, 600), ylabel="∑E/E0", titlelocation=:left, legend=false)
    title!(p_[1], "Reflection")
    title!(p_[2], "Axion Induced")
    plot!(p_[1], frequencies / 1e9, [p[2][1] for p in P]; xformatter=_ -> "", xlabel="")
    plot!(p_[2], frequencies / 1e9, [p[2][2] for p in P]; xlabel="Frequency [GHz]")

    display(p_)

    return
end


function showFields(data::Data, frequencies::Vector{Float64}, freqsplot=nothing)
    if isnothing(freqsplot)
        freqsplot = data.freqs
    end

    P = Plots.Plot[]

    for i in eachindex(data)
        append!(P, [showFields(reshape(data[i].pos, size(data[i].pos, 1)), f, freqsplot;
            eps=data.s.eps, tand=data.s.tand, obj=data.obj[i])[1]
                    for f in frequencies])
    end

    xlim1 = maximum([xlims(p[1])[2] for p in P])
    ylim1 = maximum([ylims(p[1])[2] for p in P])
    ylim2 = maximum([ylims(p[2])[2] for p in P])
    ylim3 = maximum([ylims(p[3])[2] for p in P])

    for i in eachindex(P)
        xlims!(P[i][1], (-0.002, xlim1))
        xlims!(P[i][2], (-0.002, xlim1))
        ylims!(P[i][1], (-ylim1, ylim1))
        ylims!(P[i][2], (-ylim2, ylim2))
        ylims!(P[i][3], (-0.1, ylim3))

        display(P[i])
    end

    return
end

