

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

    tag::Symbol

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

function getTag(path)
    if occursin("NM",path)
        return :nm
    elseif occursin("SA",path)
        return :sa
    elseif occursin("LS",path)
        return :ls
    else
        return Symbol()
    end
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
    
    tags::Vector{Symbol}

    s::Settings
end

import Base: getindex, length, eachindex, append!

function Base.getindex(data::Data,inds::Vector{Int64})
    return Data(
        data.pos[:,inds],
        data.dist[:,inds],
        data.boost[:,inds],
        data.ref[:,inds],
        data.freqs,
        data.obj[inds],
        data.runtime[inds],
        data.opttime[inds],
        data.optdist[inds],
        data.tags[inds],
        data.s
    )
end

function Base.getindex(data::Data,idx::Integer)
    return Data(
        data.pos[:,[idx]],
        data.dist[:,[idx]],
        data.boost[:,[idx]],
        data.ref[:,[idx]],
        data.freqs,
        data.obj[[idx]],
        data.runtime[[idx]],
        data.opttime[[idx]],
        data.optdist[[idx]],
        data.tags[[idx]],
        data.s
    )
end

Base.length(data::Data) = length(data.obj)
Base.eachindex(data::Data) = eachindex(data.obj)

function append!(database::Data,dataappend::Data)
    @assert database.s == dataappend.s

    database.pos = hcat(database.pos,dataappend.pos)
    database.dist = hcat(database.dist,dataappend.dist)
    database.boost = hcat(database.boost,dataappend.boost)
    database.ref = hcat(database.ref,dataappend.ref)
    append!(database.obj,dataappend.obj)
    append!(database.runtime,dataappend.runtime)
    append!(database.opttime,dataappend.opttime)
    append!(database.optdist,dataappend.optdist)
    append!(database.tags,dataappend.tags)
    
    return database
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

    tags = fill(getTag(path),size(data,1))

    println("$(size(data,1)-sum(idxs)) rejected out of $(size(data,1)).")

    data = Data(
        pos',
        dist',
        boost',
        ref',
        freqs,
        data[idxs,s.ndisk+1],
        T[idxs],
        data[idxs,s.ndisk+2],
        data[idxs,s.ndisk+3],
        tags,
        s
    )

    return data[unique(i -> data.obj[i], 1:length(data))]
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
        data.tags[idx],
        data.s
    )
end


function prepareDataAll1d(path,threshold=Inf64;
        f0=22.025e9,df=50e6,nf=10,ndisk=20,eps=24.0,tand=0.0,
        filterin="",filterout="",filterany=false)

    pos_ = []
    dist_ = []
    boost_ = []
    ref_ = []
    obj_ = []
    T_ = []
    opttime_ = []
    optdist_ = []
    tags_ = []

    n, r = 0, 0

    s = Settings(f0,df,nf,ndisk,eps,tand)
    pathes = []
    p = "$(f0)_$(df)_$(nf)_$(ndisk)_$(eps)_$(tand)"

    for (root, dirs, files) in walkdir(path)
        for path in joinpath.(root,files)
            if path in pathes || !occursin(".jld2",path) || !occursin(p,path)
                continue
            end

            if filterin isa AbstractArray
                if filterany
                    if !any(occursin.(filterin,path))
                        continue
                    end
                else
                    if !all(occursin.(filterin,path))
                        continue
                    end
                end
            elseif !occursin(filterin,path)
                continue
            end

            if filterout isa AbstractArray
                if filterany
                    if any(occursin.(filterout,path))
                        continue
                    end
                else
                    if all(occursin.(filterout,path))
                        continue
                    end
                end
            elseif filterout != "" && occursin(filterout,path)
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
            
            tags = fill(getTag(path),size(data,1))

            push!(pos_,pos')
            push!(dist_,dist')
            push!(boost_,boost')
            push!(ref_,ref')
            push!(obj_,data[idxs,s.ndisk+1])
            push!(T_,T[idxs])
            push!(opttime_,data[idxs,s.ndisk+2])
            push!(optdist_,data[idxs,s.ndisk+3])
            push!(tags_,tags)

            r += size(data,1)-sum(idxs); n += size(data,1)
        end
    end

    println("Data preparation: $r rejected out of $n.")
    println("Data preparation: $(n-r) accepted.")
    
    freqs = genFreqs(s.f0,s.df; n=s.nf)

    _pos =      cat(pos_...; dims=2)
    _dist =     cat(dist_...; dims=2)
    _boost =    cat(boost_...; dims=2)
    _ref =      cat(ref_...; dims=2)
    _obj =      cat(obj_...; dims=1)
    _T =        cat(T_...; dims=1)
    _opttime =  cat(opttime_...; dims=1)
    _optdist =  cat(optdist_...; dims=1)
    _tags =     cat(tags_...; dims=1)

    data = Data(
        _pos,
        _dist,
        _boost,
        _ref,
        freqs,
        _obj,
        _T,
        _opttime,
        _optdist,
        _tags,
        s
    )

    return data[unique(i -> data.obj[i], 1:length(data))]
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
    data.tags .= data.tags[sp]

    return
end




function findOutliers(data::Data,threshold::Float64;
        showdistribution::Bool=false,bymax::Bool=false)

    M = reshape(sum(data.dist,dims=2),size(data.dist,1))/size(data.dist,2)
    D = data.dist .- M

    if bymax
        d = reshape(maximum(abs.(D),dims=1),size(data.dist,2))
    else
        d = reshape(sum(abs.(D),dims=1),size(data.dist,2))
    end

    idxs = findall(x->x>threshold,d)#; push!(idxs,argmin(data.obj))

    if showdistribution
        lim = 20*round(log(10,length(idxs)))
        p1 = scatter(1:size(data.dist,1),M/1e-3;
            xlabel="Disc Index",ylabel="Distances [mm]",title="Average Position",legend=false)

        p2 = histogram(d/1e-3; ylims=(0,lim),bins=100,
            xlabel=bymax ? "Max. Distance From Average [mm]" :
                "Summed Abs. Distance From Average [mm]",
            ylabel="Counts",legend=false)
        vline!(p2,[threshold/1e-3])

        display(p1)
        display(p2)
    end

    return data[idxs]
end





# function saveData(data)
#     if is
#         df = DataFrame()

#         for i in 1:data.s.ndisk
#             df[!,"d$i"] = data.dist[i,:]
#         end
#         df[!,"obj"] = data.obj
#     end

#     sort!(df,[:obj])
#     unique!(df)

#     h5write("test.h5","data",df)

#     return
# end





# P0 = Dict{Int,Vector{Float64}}()
# B0 = Vector{Float64}([])
# fails = []

# for i in 10:100
#     try
#         data = prepareDataAll1d(getPath();f0=i*1e9+25e6,tand=6e-5,)
#         if length(data) == 0
#             error("No data found for $i GHz.")
#         end
#         b = best(data)
#         P0[i] = b.pos[:]
#         push!(B0,b.obj[1])
#     catch e
#         println("No data found for ",i," GHz.\n")
#         P0[i] = zeros(20)
#         push!(B0,0)
#         push!(fails,i)
#     end
# end

# for i in fails
#     try
#         data = prepareDataAll1d(getPath();f0=i*1e9+25e6,tand=6e-5,)
#         if length(data) == 0
#             error("No data found for $i GHz.")
#         end
#         b = best(data)
#         P0[i] = b.pos[:]
#         B0[i-9] = b.obj[1]
#     catch e
#         println(e)
#     end
# end

# @save "examples/full_20_24.0_6.0e-5.jld2" P0 B0

