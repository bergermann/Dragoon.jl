

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

    println("Data preparation: $(size(data,1)-sum(idxs)) rejected out of $(size(data,1)).")

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

    return Data(
        pos',
        dist',
        boost',
        ref',
        freqs,
        data[idxs,s.ndisk+1],
        T,
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


function prepareDataAll1d(path,threshold=Inf64)
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

    for (root, dirs, files) in walkdir(path)
        p = joinpath.(root,files)

        if length(p) != 1 || !occursin(".jld2",p[1])
            continue
        end

        @load String(p[1]) data sigx s seed T
        
        idxs = (data[:,s.ndisk+1] .<= threshold)

        r += size(data,1)-sum(idxs); n += size(data,1)

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
        push!(T_,T)
        push!(opttime_,data[idxs,s.ndisk+2])
        push!(optdist_,data[idxs,s.ndisk+3])
    end

    println("Data preparation: $r rejected out of $n.")
    
    freqs = genFreqs(s.f0,s.df; n=s.nf)

    return Data(
        cat(pos_...; dims=2),
        cat(dist_...; dims=2),
        cat(boost_...; dims=2),
        cat(ref_...; dims=2),
        freqs,
        cat(obj_...; dims=1),
        cat(T_...; dims=1),
        cat(opttime_...; dims=1),
        cat(optdist_...; dims=1),
        s
    )
end



function histogramB(data; kwargs...)
    p = histogram(data; legend=false,kwargs...)
    vline!(p,[minimum(data),maximum(data),sum(data)/length(data)])

    return p
end

function showDist(data,n=typemax(Int); kwargs...)
    n = min(n,size(data.dist,2))

    p = scatter(1:size(data.dist,1),data.dist[:,1:n];
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
    hline!(p1,[-threshold])
    display(p1)

    p2 = histogram2d(data.optdist,-data.obj; legend=false,
        xlabel="travel distance [m]",ylabel="objective value")
    hline!(p2,[-threshold])
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