


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



