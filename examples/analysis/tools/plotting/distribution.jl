




function showDist(data,n=typemax(Int); kwargs...)
    n = min(n,size(data.dist,2))

    p = scatter(1:size(data.dist,1),data.dist[:,1:n]*1e3;
        legend=false,c=:blue,markersize=2,kwargs...)

    return p
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

    ridge(Hdx,Hdy; normalize=true,yscale=1.5,xlabel="Δd [mm]",ylabel="Disk index i")
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
    
    ridge(Hpx,Hpy; normalize=true,yscale=1.5,xlabel="Δp [mm]",ylabel="Disk index i")
    display(vline!([0]; c=:black,linestyle=:dash,linewidth=0.5))

    mp = reshape(mean(dp; dims=1),(size(dd,2)));
    display(histogram(mp/1e-3; xlabel="μ(Δp) [mm]",bins=-0.1:0.001:0.1,label=""))

    return
end

