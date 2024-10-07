

function plotRescale(booster,p0,cf,bw,df,s,n1,n2; ylims=(-1,39))
    move(booster,p0; additive=false)
    f = genFreqs(cf,3*bw; n=1000)

    p1 = plot(f/1e9,getBoost1d(booster,f)/1e3;
        xlabel="Frequency [GHz]",ylabel="Boostfactor β² × 10³",c=:red,label="original",lw=2,ylims=ylims)
    p2 = plot(f/1e9,abs.(getRef1d(booster,f));
        xlabel="Frequency [GHz]",ylabel="Reflectivity |R|",c=:red,label="original",lw=2)

    for i in n1:n2
        if i == 0
            continue
        end
        
        scale = cf/(cf+i*df)
        dd = (scale-1)*pos2dist(p0; disk_thickness=booster.thickness)

        move(booster,dist2pos(pos2dist(p0)+dd*s); additive=false)
        
        f = genFreqs(cf+df*i,3*bw; n=1000)
        plot!(p1,f/1e9,getBoost1d(booster,f)/1e3; c=:blue,label=(i == n1 ? "rescaled" : ""))
        plot!(p2,f/1e9,abs.(getRef1d(booster,f)); c=:blue,label=(i == n1 ? "rescaled" : ""))
    end

    display(p1)
    display(p2)

    return p1, p2
end