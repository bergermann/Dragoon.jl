

mutable struct Settings
    f0::Float64
    df::Float64
    nf::Int
    ndisk::Int
    eps::Float64
    tand::Float64
end

struct Entry
    pos::Vector{Float64}
    dist::Vector{Float64}

    boost::Vector{Float64}
    ref::Vector{ComplexF64}

    obj::Float64

    runtime::Float64
    opttime::Float64
    optdist::Float64

    s::Settings
end

function getPath(sigx,Nsig,f0,df,nf,ndisk,eps,tand,algorithm="",time="")
    p = joinpath(
        "examples",
        "analysis",
        "optimization data",
        "$(sigx)_$(Nsig)_$(f0)_$(df)_$(nf)_$(ndisk)_$(eps)_$(tand)",
        algorithm,
        isempty(time) ? "" : time*".jld2"
    )

    @assert ispath(p) "No directory or file exists at $p"

    return p
end

getPath() = joinpath("examples","analysis","optimization data")



mutable struct Data
    pos::Matrix{Float64}
    dist::Matrix{Float64}

    boost::Matrix{Float64}
    ref::Matrix{ComplexF64}

    freqs::Vector{Float64}

    obj::Vector{Float64}

    runtime::Vector{Float64}
    opttime::Vector{Float64}
    optdist::Vector{Float64}

    s::Settings
end

function prepareData1d(path,threshold=Inf64)
    @load String(path) data sigx s seed T
    
    idxs = (data[:,s.ndisk+1] .<= threshold)

    print("Data preparation: ")

    idxs = (data[:,s.ndisk+1] .<= threshold)
    idxs_ = all(data[:,1:s.ndisk-1] .<= data[:,2:s.ndisk],dims=2)
    idxs .*= reshape(idxs_,(length(idxs_),))

    pos = data[idxs,1:s.ndisk]
    dist = similar(pos)

    for i in axes(pos,1)
        dist[i,:] = pos2dist(pos[i,:])
    end

    freqs = genFreqs(s.f0,s.df; n=s.nf)

    boost = zeros(Float64,size(dist,1),s.nf)
    ref = zeros(ComplexF64,size(dist,1),s.nf)

    for i in axes(dist,1)
        boost[i,:] = boost1d(dist[i,:],freqs; eps=s.eps,tand=s.tand)
        ref[i,:] = ref1d(dist[i,:],freqs; eps=s.eps,tand=s.tand)
    end

    println("$(size(data,1)-sum(idxs)) rejected out of $(size(data,1)).")

    return Data(
        pos',
        dist',
        boost',
        ref',
        freqs,
        data[idxs,s.ndisk+1],
        T[idxs],
        data[idxs,s.ndisk+2],
        data[idxs,s.ndisk+3],
        s
    )
end

function best(data::Data)
    idx = argmin(data.obj)

    return Entry(
        data.pos[:,idx],
        data.dist[:,idx],
        data.boost[:,idx],
        data.ref[:,idx],
        data.obj[idx],
        data.runtime[idx],
        data.opttime[idx],
        data.optdist[idx],
        data.s
    )
end


function prepareDataAll1d(path,threshold=Inf64; f0=22.025e9,df=50e6,nf=10,ndisk=20,eps=24.0,tand=0.0)
    pos_ = []
    dist_ = []
    boost_ = []
    ref_ = []
    obj_ = []
    T_ = []
    opttime_ = []
    optdist_ = []

    n, r = 0, 0

    s = 0
    pathes = []
    p = "$(f0)_$(df)_$(nf)_$(ndisk)_$(eps)_$(tand)"

    for (root, dirs, files) in walkdir(path)
        for path in joinpath.(root,files)
            if path in pathes || !occursin(".jld2",path) || !occursin(p,path)
                continue
            end

            push!(pathes,path)

            println("opening ",String(path))
            @load String(path) data sigx s seed T
            
            idxs = (data[:,s.ndisk+1] .<= threshold)
            idxs_ = all(data[:,1:s.ndisk-1] .<= data[:,2:s.ndisk],dims=2)
            idxs .*= reshape(idxs_,(length(idxs_),))

            pos = data[idxs,1:s.ndisk]
            dist = similar(pos)

            for i in axes(pos,1)
                dist[i,:] = pos2dist(pos[i,:])
            end

            freqs = genFreqs(s.f0,s.df; n=s.nf)

            boost = zeros(Float64,size(dist,1),s.nf)
            ref = zeros(ComplexF64,size(dist,1),s.nf)

            for i in axes(dist,1)
                boost[i,:] = boost1d(dist[i,:],freqs; eps=s.eps,tand=s.tand)
                ref[i,:] = ref1d(dist[i,:],freqs; eps=s.eps,tand=s.tand)
            end

            push!(pos_,pos')
            push!(dist_,dist')
            push!(boost_,boost')
            push!(ref_,ref')
            push!(obj_,data[idxs,s.ndisk+1])
            push!(T_,T[idxs])
            push!(opttime_,data[idxs,s.ndisk+2])
            push!(optdist_,data[idxs,s.ndisk+3])

            r += size(data,1)-sum(idxs); n += size(data,1)
        end
    end

    println("Data preparation: $r rejected out of $n.")
    
    freqs = genFreqs(s.f0,s.df; n=s.nf)

    _pos =      cat(pos_...; dims=2)
    _dist =     cat(dist_...; dims=2)
    _boost =    cat(boost_...; dims=2)
    _ref =      cat(ref_...; dims=2)
    _obj =      cat(obj_...; dims=1)
    _T =        cat(T_...; dims=1)
    _opttime =  cat(opttime_...; dims=1)
    _optdist =  cat(optdist_...; dims=1)

    return Data(
        _pos,
        _dist,
        _boost,
        _ref,
        freqs,
        _obj,
        _T,
        _opttime,
        _optdist,
        s
    )
end




function sortData!(data::Data)
    sp = sortperm(data.obj)

    data.pos .= data.pos[:,sp]
    data.dist .= data.dist[:,sp]
    data.boost .= data.boost[:,sp]
    data.ref .= data.ref[:,sp]

    data.obj .= data.obj[sp]
    data.runtime .= data.runtime[sp]
    data.opttime .= data.opttime[sp]
    data.optdist .= data.optdist[sp]

    return
end




function histogramB(data; kwargs...)
    p = histogram(data; legend=false,kwargs...)
    vline!(p,[minimum(data),maximum(data),sum(data)/length(data)])

    return p
end

function showDist(data,n=typemax(Int); kwargs...)
    n = min(n,size(data.dist,2))

    p = scatter(1:size(data.dist,1),data.dist[:,1:n]*1e3;
        legend=false,c=:blue,markersize=2,kwargs...)

    return p
end

function showClusters(data,k,showmax=typemax(Int))
    showmax = min(showmax,size(data.dist,2))
    ndisk = size(data.dist,1)

    c = kmeans(data.dist,k)

    colors = palette([:red,:blue],k)
    ylim = [0.98*minimum(data.dist),1.02*maximum(data.dist)]

    display(scatter(1:ndisk,data.dist[:,1:showmax];
        c=colors[c.assignments'],legend=false,ylims=ylim,title="all",markersize=2))
    display(scatter(1:ndisk,c.centers;
        c=colors,legend=false,ylims=ylim,title="centers - dist",markersize=2))

    for i in 1:k
        idxs = (1:size(data.dist,2))[c.assignments .== i]
        idxs = idxs[1:min(length(idxs),showmax)]

        display(scatter(1:ndisk,data.dist[:,idxs];
            c=colors[i],legend=false,ylims=ylim,title="assignments $i - dist",markersize=2))
    end

    pbc = plot(; legend=false,title="centers - boost",seriestype=:line)
    for i in 1:k
        boost = boost1d(c.centers[:,i],data.freqs; eps=data.s.eps,tand=data.s.tand)
        plot!(pbc,data.freqs/1e9,boost; c=colors[i])
    end
    display(pbc)

    for i in 1:k
        idxs = (1:size(data.dist,2))[c.assignments .== i]
        idxs = idxs[1:min(length(idxs),showmax)]

        display(plot(data.freqs/1e9,data.boost[:,idxs];
            c=colors[i],legend=false,title="assignments $i - boost"))
    end

    prc = plot(; legend=false,title="centers - ref",seriestype=:line)
    for i in 1:k
        ref = ref1d(c.centers[:,i],data.freqs; eps=data.s.eps,tand=data.s.tand)
        plot!(prc,data.freqs/1e9,real.(ref); c=colors[i],linestyle=:solid)
        plot!(prc,data.freqs/1e9,imag.(ref); c=colors[i],linestyle=:dash)
    end
    display(prc)

    for i in 1:k
        idxs = (1:size(data.dist,2))[c.assignments .== i]
        idxs = idxs[1:min(length(idxs),showmax)]

        pra = plot(data.freqs/1e9,real.(data.ref[:,idxs]);
            c=colors[i],linestyle=:solid,legend=false,title="assignments $i - ref")
        plot!(pra,data.freqs/1e9,imag.(data.ref[:,idxs]);
            c=colors[i],linestyle=:dash)

        display(pra)
    end

    return c
end

function showQuality(data,threshold=nothing)
    display(histogramB(data.obj; title="objective values",
        xlabel="objective value"))
    display(histogramB(data.opttime; title="total travel time",
        xlabel="travel time [s]"))
    display(histogramB(data.optdist; title="total travel distance",
        xlabel="trave distance [m]"))

    p1 = histogram2d(data.opttime,-data.obj; legend=false,
        xlabel="travel time [s]",ylabel="objective value")
    if !isnothing(threshold)
        hline!(p1,[-threshold])
    end
    display(p1)

    p2 = histogram2d(data.optdist,-data.obj; legend=false,
        xlabel="travel distance [m]",ylabel="objective value")
    if !isnothing(threshold)
        hline!(p2,[-threshold])
    end
    display(p2)

    p3 = histogram2d(data.optdist,data.opttime; legend=false,
        xlabel="travel distance [m]",ylabel="travel time [m]")
    display(p3)

    if !isnothing(threshold)
        sr = round(100*sum((data.obj .<= threshold))/length(data.obj); sigdigits=4)
        println("success rate: $sr")
    end

    return
end




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




function showDistribution(data)
    dd = data.dist .- d0;
    for i in 1:20
        display(histogram(dd[i,:]/1e-3; xlabel="Δd_$i [mm]",xlim=[-1,1],ylim=[-500,40_000],
            bins=-1:0.01:1))
    end

    md = reshape(mean(dd; dims=1),(size(dd,2)));
    display(histogram(md/1e-3; xlabel="μ(Δd) [mm]",bins=-0.1:0.001:0.1))


    dp = data.pos .- p0;
    for i in 1:20
        display(histogram(dp[i,:]/1e-3; xlabel="Δp_$i [mm]",xlim=[-1,1],ylim=[-500,40_000],
            bins=-1:0.01:1))
    end

    mp = reshape(mean(dp; dims=1),(size(dd,2)));
    display(histogram(mp/1e-3; xlabel="μ(Δp) [mm]",bins=-0.1:0.001:0.1))

    return
end






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

    labels = reshape(collect(string.(round.(Int,q*100)).*"%"),(1,ndiv+1))

    p1 = plot(sigxs/1e-6,Q; xlabel="σ_Δp [μm]",ylabel="Boostfactor β²",label=labels)

    display(p1)

    return
end




