export dragoon, rescale

function dragoon(booster::Booster,hist::Vector{State},bandwidth::Float64,overlap::Real,
        objective::Callback,unstuckinator::Callback;
        fmin::Float64=10e9,fmax::Float64=100e9,nfreqs::Int=10,niters::Int=1_000,
        scalerange::NTuple{2,Float64}=(1,0,1.3),scalesteps::Int=100,
        preoptimize::Bool=true,reset::Bool=false,reverse::Bool=false)

    @assert fmin < fmax "Maximum frequency needs to be higher than minimum frequency."
    @assert nfreqs > 1 "Need at least 2 frequency points nfreqs."
    @assert bandwidth > overlap "Overlap needs to be smaller than bandwidth."
    @assert scalesteps > 1
    @assert scalerange[2] > scalerange[1] > 0

    if reverse
        freqs = collect(range(fmax-bandwidth,fmax,nfreqs))
    else
        freqs = collect(range(fmin,fmin+bandwidth,nfreqs))
    end

    if preoptimize
        if reset
            initdist = findpeak1d(22.025e9,booster.ndisk)

            move(booster,dist2pos(ones(booster.ndisk)*initdist); additive=false)
        end

        trace = nelderMead(booster,hist,freqs,
                    1.,1+2/booster.ndisk,0.75-1/2booster.ndisk,1-1/booster.ndisk,1e-12,1e-12,
                    objective,
                    InitSimplexRegular(1e-4),
                    DefaultSimplexSampler,
                    UnstuckDont;
                    maxiter=Int(1e5),
                    showtrace=false,
                    unstuckisiter=true,)

        obj_pre = updateHist!(booster,hist,freqs,objective; force=true)

        println("Preoptimization complete with objective value $(round(obj_pre; digits=2))")
    end

    i = 1; i_ = ceil(Int,(fmax-fmin)/(bandwidth-overlap))+1; t1 = copy(booster.timestamp)

    Obj = Float64[]; Pos = Vector{Float64}[]; Freqs = Vector{Float64}[]
    S = Float64[];

    cont = true

    while cont
        trace = nelderMead(booster,hist,freqs,
                    1.,1+2/booster.ndisk,0.75-1/2booster.ndisk,1-1/booster.ndisk,1e-12,1e-12,
                    objective,
                    InitSimplexRegular(1e-5),
                    DefaultSimplexSampler,
                    unstuckinator;
                    maxiter=niters,
                    showtrace=false,
                    showevery=100,
                    unstuckisiter=true,
                    resettimer=false)
        
        push!(Obj,updateHist!(booster,hist,freqs,objective))
        push!(Pos,copy(booster.pos),)
        push!(Freqs,copy(freqs))
        
        if reverse
            cont = freqs[end] > fmin
            freqs = collect(range(fmax-bandwidth*(i+1)+overlap*i,fmax-(bandwidth-overlap)*i,nfreqs))

            cf = (freqs[1]+freqs[end])/2
            scale = (cf+(bandwidth-overlap))/cf
        else
            cont = freqs[1] < fmax
            freqs = collect(range(fmin+(bandwidth-overlap)*i,fmin+bandwidth*(i+1)-overlap*i,nfreqs))

            cf = (freqs[1]+freqs[end])/2
            scale = (cf-(bandwidth-overlap))/cf
        end

        s = rescale(booster,hist,freqs,objective,scale,scalerange,scalesteps)

        push!(S,s)

        println("finished iteration $i/$i_")

        i += 1
    end

    t2 = copy(booster.timestamp)

    println("Elapsed movement time: ",canonicalize(t2-t1))

    return Obj, Pos, Freqs, S
end


function rescale(booster::Booster,hist::Vector{State},freqs::Array{Float64},objective::Callback,
        scale::Float64,scalerange::Tuple{Float64,Float64},scalesteps::Int)

    p_ = copy(booster.pos)
    dd = (scale-1)*pos2dist(booster.pos; disk_thickness=booster.thickness)

    b_ = 0; B = zeros(Float64,scalesteps+1)
    i_ = 0

    for i in 0:scalesteps
        move(booster,dist2pos(pos2dist(p_)+dd*lerp(scalerange,i/scalesteps)); additive=false)
        updateHist!(booster,hist,freqs,objective)
        
        b = sum(getBoost1d(booster,freqs))
        B[i+1] = b

        if b > b_
            b_ = b
            i_ = i
        end
    end

    # display(plot(lerp.(scalerange[1],scalerange[2],(0:scalesteps)./scalesteps),B,
    #     title="rescaling, $(lerp(scalerange,i_/scalesteps))"))

    move(booster,dist2pos(pos2dist(p_)+dd*lerp(scalerange,i_/scalesteps)); additive=false)

    return lerp(scalerange,i_/scalesteps)
end