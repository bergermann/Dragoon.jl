
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

    s
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

    pbc = plot(; legend=false,ylims=ylim,title="centers - boost",markersize=2)
    for i in 1:k
        scatter!(pbc,1:ndisk,
            boost1d(c.centers[:,i],data.freqs; eps=data.s.eps,tand=data.s.tand);
            c=colors[i])
    end
    display(pbc)

    prc = plot(; legend=false,ylims=ylim,title="centers - ref")
    for i in 1:k
        ref = ref1d(c.centers[:,i],data.freqs; eps=data.s.eps,tand=data.s.tand)
        plot!(prc,1:ndisk,real.(ref); c=colors[i],linestyle=:solid)
        plot!(prc,1:ndisk,imag.(ref); c=colors[i],linestyle=:dash)
    end
    # display(prc)

    for i in 1:k
        idxs = (1:size(data.dist,2))[c.assignments .== i]
        idxs = idxs[1:min(length(idxs),showmax)]

        display(plot(data.freqs,data.boost[:,idxs];
            c=colors[i],legend=false,title="assignments $i - boost"))
    end

    return
end