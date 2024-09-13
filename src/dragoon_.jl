export dragoon

function dragoon(booster::Booster,hist::Vector{State},bandwidth::Float64,overlap::Real,
        objective::Callback,unstuckinator::Callback;
        fmin::Float64=10e9,fmax::Float64=100e9,nfreqs::Int=10,
        scalerange::NTuple{2,Int}=(0.7,1.3),scalesteps::Int=100)

    @assert fmax > fmin "Maximum frequency needs to be higher than minimum frequency."
    @assert nfreqs > 1 "Need at least 2 frequency points nfreqs."
    @assert bandwidth > overlap "Overlap needs to be smaller than bandwidth."
    @assert scalesteps > 1
    @assert scalerange[2] > scalerange[1] > 0

    freqs = range(fmin,fmin+bandwidth,nfreqs)

    i = 1

    while freqs[1] < fmax
        trace = nelderMead(booster,hist,freqs,
                    1.,1+2/booster.ndisk,0.75-1/2booster.ndisk,1-1/booster.ndisk,1e-9,
                    objective,
                    InitSimplexRegular(1e-5),
                    DefaultSimplexSampler,
                    unstuckinator;
                    maxiter=Int(1e3),
                    showtrace=true,
                    showevery=100,
                    unstuckisiter=true,)

        display(plot(freqs/1e9,getBoost1d(booster,freqs)))
        
        freqs = range(fmin+(bandwidth-overlap)*i,fmin+bandwidth*(i+1)-overlap*i,nfreqs)
        i += 1

        scale = freqs[1]/(freqs[1]-(bandwidth-overlap))

        rescale(booster,hist,freqs,obj,scale,scalerange,scalesteps)

    end

    return
end

function rescale(booster::Booster,hist::Vector{State},freqs::Array{Float64},obj::Callback,
        scale::Float64,scalerange::Tuple{Float64,Float64},scalesteps::Int)

    p0 = copy(booster.pos)
    dd = (scale-1)*pos2dist(booster.pos; disk_thickness=booster.τ)

    b_ = 0; B = zeros(Float64,scalesteps+1)
    i_ = 0

    freqsplot = range(2*freqs[1]-freqs[end],2*freqs[2]-freqs[1],100)

    p = plot(freqsplot/1e9,getBoost1d(booster,freqsplot))

    for i in 0:scalesteps
        move(booster,p0+dd*lerp(scalerange,i/scalesteps))
        updateHist!(booster,hist,freqs,obj)
        b = sum(getBoost1d(booster,freqs))

        if b > b_
            b_ = b
            i_ = i
        end

        println("step:  ",i,"/",scalesteps)
        println("scale: ",lerp(scalerange,i/scalesteps))
    end

    move(booster,p0+dd*lerp(scalerange,i_/scalesteps))

    return
end