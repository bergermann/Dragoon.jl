




function showDist(data,n=1000; ndiv=10)
    @assert iseven(ndiv) "Divisions ndiv need to be even."

    n = min(n,size(data.dist,2))

    p1 = scatter(1:size(data.dist,1),data.dist[:,1:n]*1e3;
        legend=false,c=:blue,markersize=2,xlabel="Disc Index",ylabel="Distance [m]")

    q = collect(0:1/ndiv:1)
    Q = zeros(Float64,data.s.ndisk,ndiv+1)

    for i in 1:data.s.ndisk
        Q[i,:] = quantile(data.dist[i,:],q)/1e-3
    end

    display(Q)

    m = mean(data.dist,dims=2)/1e-3

    p2 = plot(; xlabel="Disc Index",ylabel="Disc Distance [mm]",
        title="Distance Quantiles of Optimized States")

    for i in 1:div(length(q) - 1, 2)
        plot!(p2,1:data.s.ndisk,Q[:,i]; linewidth=0,
            fillrange=Q[:,end-i+1],fillalpha=0.3,fillcolor=1,label=""
        )
    end

    scatter!(p2,1:data.s.ndisk,m; markersize=1,label="mean")

    # display(p1)
    display(p2)

    return 
end


function showDistribution(data,d0)
    dd = data.dist .- d0
    Hdx = []; Hdy = []

    for i in 1:data.s.ndisk
        f = fit(Histogram,dd[i,:]/1e-3,-1:0.01:1)
        e = collect(f.edges[1])
        push!(Hdx,(e[2:end]+e[1:end-1])/2)
        push!(Hdy,f.weights)
    end

    ridge(Hdx,Hdy; normalize=true,yscale=1.5,xlabel="Δd [mm]",ylabel="Disk index i",
        title="Deviation from Initial Config in Distance Space")
    display(vline!([0]; c=:black,linestyle=:dash,linewidth=0.5))

    md = reshape(mean(dd; dims=1),(size(dd,2)));
    display(histogram(md/1e-3; xlabel="μ(Δd) [mm]",bins=-0.1:0.001:0.1,label=""))

    dp = data.pos .- p0;
    Hpx = []; Hpy = []

    for i in 1:data.s.ndisk
        f = fit(Histogram,dp[i,:]/1e-3,-1:0.01:1)
        e = collect(f.edges[1])
        push!(Hpx,(e[2:end]+e[1:end-1])/2)
        push!(Hpy,f.weights)
    end
    
    ridge(Hpx,Hpy; normalize=true,yscale=1.5,xlabel="Δp [mm]",ylabel="Disk index i",
        title="Deviation from Initial Config in Position Space")
    display(vline!([0]; c=:black,linestyle=:dash,linewidth=0.5))

    mp = reshape(mean(dp; dims=1),(size(dd,2)));
    display(histogram(mp/1e-3; xlabel="μ(Δp) [mm]",bins=-0.1:0.001:0.1,label=""))

    return
end

