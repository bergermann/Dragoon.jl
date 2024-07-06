
using Plots, ColorSchemes, Clustering, BoostFractor, Dragoon, Statistics, Peaks
using JLD2, HDF5, DataFrames

include("tools.jl");

# path = getPath(1e-3,100_000,22.025e9,50e6,10,20,24.0,0.0);
path = getPath("rand",100_000,22.025e9,50e6,10,20,24.0,0.0,"NM2","2024_05_29-15_05_16");

# data = prepareDataAll1d(path,-10_000);
data = prepareDataAll1d(getPath(),-14_000);
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

showQuality(data)

showDist(data,1000; xlabel="disc index", ylabel="d_i [mm]")
hline!([initdist]*1e3)

# plot(data.freqs/1e9,data.boost; xlabel="frequency [GHz]",ylabel="boost factor β^2")

# c = showClusters(data,100,200);
# scanSingleDiscs(best(data).dist,22.025e9,50e6,50,250,100)

showDistribution(data)



function showFields(pos::Vector{Float64},frequency::Float64;
        R::Float64=0.1,eps::Float64=24.,tand::Float64=0.,thickness::Float64=1e-3)
    
    ndisk = length(pos)

    # eps_ = ComplexF64[!iseven(i) ? 1 : 24 for i in 1:2*ndisk+1]
    # prepend!(eps_,NaN)

    epsilon = 24
    eps_ = Array{Complex{Float64}}([NaN, 1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1,epsilon,1])
    

    # d = pos2dist(pos; disk_thickness=thickness)
    # distance = Float64[0]
    # for i in 1:ndisk; append!(distance,d[i],thickness); end
    # append!(distance,0)
    distance = [0, 1.00334, 1.0,
                    6.94754, 1.0,
                    7.17660, 1.0,
                    7.22788, 1.0,
                    7.19717, 1.0,
                    7.23776, 1.0,
                    7.07746, 1.0,
                    7.57173, 1.0,
                    7.08019, 1.0,
                    7.24657, 1.0,
                    7.21708, 1.0,
                    7.18317, 1.0,
                    7.13025, 1.0,
                    7.21980, 1.0,
                    7.45585, 1.0,
                    7.39873, 1.0,
                    7.15403, 1.0,
                    7.14252, 1.0,
                    6.83105, 1.0,
                    7.42282, 1.0,
                    0.0]*1e-3

    coords = SeedCoordinateSystem(X = [1e-9], Y = [1e-9])
    
    sbdry = SeedSetupBoundaries(coords; diskno=ndisk,distance=distance,epsilon=eps_)
    
    modes = SeedModes(coords, ThreeDim=false, Mmax=1, Lmax=0, diskR=R)
    m_reflect = Float64[1.0]

    boost, refl = transformer(sbdry,coords,modes;
        reflect=m_reflect,prop=propagator1D,diskR=R,f=frequency)
    Eout = transpose([boost refl])

    full_fields = BoostFractor.transformer_trace_back(refl,m_reflect, sbdry, coords,modes;
        prop=propagator1D,f=frequency)

    # return full_fields[:,:,1]

    plot_1d_field_pattern(-autorotate(full_fields[:,:,1]), sbdry, frequency)

    
    # full_fields = BoostFractor.transformer_trace_back(EoutModes0[1,:,freq_idx],
    #     zeros(length(m_reflect)), sbdry, coords,modes;
    #     prop=propagator1D,f=frequencies[freq_idx], inlcudes_axion=true)
    # plot_1d_field_pattern(-(full_fields[:,:,1]), sbdry, frequencies[freq_idx],
    #     add_ea=true, overallphase=exp(-1im*pi/2*0.95))
end


function plot_1d_field_pattern(full_solution_regions, bdry::SetupBoundaries, f; fill=false,
        add_ea=false, overallphase=1)
    
    ztot = 0 # Each region needs to know where it starts, so iteratively add up the lengths of regions
    Nregions = length(bdry.eps) 
    c = 299792458.

    p = plot()
    for s in 1:Nregions
        # Propagation constant in that region
        kreg = 2pi/c*f
        kreg *= sqrt(bdry.eps[s])
        
        # Define the color of that region according to mirror / disk / free space / ...
        fillcolor = nothing
        if abs.(bdry.eps[s]) > 100 || bdry.eps[s] == NaN
            fillcolor = "darkorange"
        elseif abs.(bdry.eps[s]) != 1
            fillcolor = "lightgray"
        end
        
        # Plot that region        
        plotRegion1d!(p,full_solution_regions[s,1].*overallphase,
                        full_solution_regions[s,2].*overallphase,
                        ztot, ztot+bdry.distance[s],
                        kreg,
                        Ea=(add_ea ? (1/bdry.eps[s]).*overallphase : 0),
                        maxE=2.2*maximum(abs.(full_solution_regions[:,:])),
                        bgcolor=fillcolor, extraspace=(s == Nregions),fill=fill,)
        
        ztot += bdry.distance[s]
    end

    return p
end

function plotRegion1d!(p,R,L,z0,z1,k; 
        bgcolor=nothing,extraspace=false,fill=false,Ea=0,maxE=10)    

    # Construct the relative coordinate system for that region
    z1 += 1e-9
    maximum = (z1+(extraspace ? 10e-3 : 0))
    z = vcat(z0:2e-4:maximum, maximum)
    dz = z .- z0
    dz2 = .-(z .- z1)
    
    # Calculate the functional solution for the region
    Rf = L*exp.(+1im.*k.*dz2)
    Lf = R*exp.(-1im.*k.*dz2)
    
    #Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Mark the Region as Disk / Mirror / Air / etc.
    # if !isnothing(bgcolor)
    #     fill_between([z0 == 0 ? -0.0025 : z0,z1], -maxE, maxE, color=bgcolor, linewidth=0)
    # end
    
    # Plot the real and imaginary part of the solution
    plot!(p,z, real.(Rf.+Lf.+ Ea), c=:blue, label=(extraspace ? "Re(E)" : ""))
    plot!(p,z, imag.(Rf.+Lf.+ Ea), c=:red, label=(extraspace ? "Im(E)" : ""))

    return
end

function autorotate(full_solution_regions)
    ang = angle.(full_solution_regions[2,1].+full_solution_regions[2,2])
    #sgn = real.((full_solution_regions[2,1].+full_solution_regions[2,2]).*exp(-1im*ang)) .> 0
    return full_solution_regions.*exp(-1im*ang)
end