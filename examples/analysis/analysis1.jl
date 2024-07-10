
using Plots, ColorSchemes, Clustering, BoostFractor, Dragoon, Statistics, Peaks
using JLD2, HDF5, DataFrames

include("tools.jl");

# path = getPath(1e-3,100_000,22.025e9,50e6,10,20,24.0,0.0);
path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0,"NM2","2024_05_29-15_05_16");

# data = prepareDataAll1d(path,-10_000);
# data = prepareDataAll1d(getPath(),-14_000);
data = prepareData1d(path,-14_000);



function saveData(data)
    if is
        df = DataFrame()

        for i in 1:data.s.ndisk
            df[!,"d$i"] = data.dist[i,:]
        end
        df[!,"obj"] = data.obj
    end

    sort!(df,[:obj])
    unique!(df)

    h5write("test.h5","data",df)

    return
end













initdist = findpeak1d(data.s.f0,20)
d0 = initdist*ones(data.s.ndisk); p0 = dist2pos(d0);
freqsplot = genFreqs(22.025e9,150e6; n=1000);

sortData!(data)
b = best(data)

# showQuality(data)

# showDist(data,1000; xlabel="disc index", ylabel="d_i [mm]")
# hline!([initdist]*1e3)

# # plot(data.freqs/1e9,data.boost; xlabel="frequency [GHz]",ylabel="boost factor β^2")

# # c = showClusters(data,100,200);
# # scanSingleDiscs(best(data).dist,22.025e9,50e6,50,250,100)

# showDistribution(data)



function showFields(pos::Vector{Float64},frequency::Float64;
        R::Float64=0.1,eps::Float64=24.,tand::Float64=0.,thickness::Float64=1e-3)
    
    ndisk = 20#length(pos)

    eps_ = ComplexF64[!iseven(i) ? 1 : eps for i in 1:2*ndisk+1]
    prepend!(eps_,NaN)

    d = pos2dist(pos; disk_thickness=thickness)
    distance = Float64[0]
    for i in 1:ndisk; append!(distance,d[i],thickness); end
    append!(distance,0)

    coords = SeedCoordinateSystem(X = [1e-9], Y = [1e-9])
    
    sbdry = SeedSetupBoundaries(coords; diskno=ndisk,distance=distance,epsilon=eps_)
    
    modes = SeedModes(coords, ThreeDim=false, Mmax=1, Lmax=0, diskR=R)
    m_reflect = Float64[1.0]

    boost, refl = transformer(sbdry,coords,modes;
        reflect=m_reflect,prop=propagator1D,diskR=R,f=frequency)
    Eout = transpose([boost refl])

    p = plot(; layout=(2,1),size=(800,600),xlabel="z Position [m]",ylabel="E/E0")
    title!(p[1],"Reflection")
    title!(p[2],"Axion Induced")
        # title="$(round(frequency/1e9,sigdigits=5)) MHz")

    drawBdry!(p[1],sbdry)
    drawBdry!(p[2],sbdry)

    full_fields_r = BoostFractor.transformer_trace_back(Eout[2,:],m_reflect, sbdry, coords,modes;
        prop=propagator1D,f=frequency)
    # plot_1d_field_pattern!(p[1],-autorotate(full_fields_r[:,:,1]), sbdry, frequency)
    plot_1d_field_pattern!(p[1],full_fields_r[:,:,1],sbdry,frequency)

    
    full_fields_a = BoostFractor.transformer_trace_back(Eout[1,:],Float64[0],sbdry,coords,modes;
        prop=propagator1D,f=frequency, inlcudes_axion=true)
    # plot_1d_field_pattern(-(full_fields_a[:,:,1]), sbdry, frequency,
    #     add_ea=true, overallphase=exp(-1im*pi/2*0.95))
    plot_1d_field_pattern!(p[2],full_fields_a[:,:,1], sbdry, frequency,
        add_ea=true, overallphase=exp(-1im*pi/2*0.95))

    return p
end


function plot_1d_field_pattern!(p,full_solution_regions, bdry::SetupBoundaries, f;
        add_ea=false,overallphase=1)
    
    ztot = 0
    Nregions = length(bdry.eps) 
    c0 = 299792458.

    maxE = maximum(abs.(full_solution_regions[:,:]))
    ylims!(p,Tuple([-1,1]*max(ylims(p)[2],2.2*maxE)))

    for s in 1:Nregions
        kreg = 2pi/c0*f*sqrt(bdry.eps[s])

        plotRegion1d!(p,full_solution_regions[s,1].*overallphase,
                        full_solution_regions[s,2].*overallphase,
                        ztot,ztot+bdry.distance[s],kreg;
                        extraspace=(s==Nregions),
                        Ea=(add_ea ? (1/bdry.eps[s]).*overallphase : 0))
        
        ztot += bdry.distance[s]
    end

    return
end

function plotRegion1d!(p,R_,L_,z0,z1,k; extraspace=false,Ea=0)    
    z1 += 1e-9
    maximum = (z1+(extraspace ? 10e-3 : 0))
    z = vcat(z0:2e-4:maximum, maximum)
    dz = z .- z0
    dz2 = .-(z .- z1)
    
    Rf = L_*exp.(+1im.*k.*dz2)
    Lf = R_*exp.(-1im.*k.*dz2)

    E = Rf.+Lf.+ Ea
    
    # plot!(p,z, real.(E); c=c,label=(extraspace ? label : ""))
    plot!(p,z, abs.(E), c=:black,linestyle=:dash,linewidth=0.5,
        label=(extraspace ? "|E|" : ""))
    plot!(p,z, real.(E), c=:blue,label=(extraspace ? "Re(E)" : ""))
    plot!(p,z, imag.(E), c=:red,label=(extraspace ? "Im(E)" : ""))
    # plot!(p,z, angle.(E), c=:black, label=(extraspace ? "∠(E)" : ""))

    return
end

function drawBdry!(p,bdry)
    z0 = 0

    for i in eachindex(bdry.eps)
        z1 = z0+bdry.distance[i]
        if abs.(bdry.eps[i]) > 100 || bdry.eps[i] == NaN
            vspan!(p,[z0 == 0 ? -0.0025 : z0,z1],c=:darkorange,label="")
        elseif abs.(bdry.eps[i]) != 1
            vspan!(p,[z0 == 0 ? -0.0025 : z0,z1],c=:lightgray,label="")
        end
    
        z0 = z1
    end

    return
end

function autorotate(full_solution_regions)
    ang = angle.(full_solution_regions[2,1].+full_solution_regions[2,2])
    #sgn = real.((full_solution_regions[2,1].+full_solution_regions[2,2]).*exp(-1im*ang)) .> 0
    return full_solution_regions.*exp(-1im*ang)
end