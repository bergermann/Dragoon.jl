export dragoon, rescale

function dragoon(booster::Booster,hist::Vector{State},bandwidth::Float64,overlap::Real,
        objective::Callback,unstuckinator::Callback;
        fmin::Float64=10e9,fmax::Float64=100e9,nfreqs::Int=10,
        scalerange::NTuple{2,Float64}=(1,0,1.3),scalesteps::Int=100)

    @assert fmax > fmin "Maximum frequency needs to be higher than minimum frequency."
    @assert nfreqs > 1 "Need at least 2 frequency points nfreqs."
    @assert bandwidth > overlap "Overlap needs to be smaller than bandwidth."
    @assert scalesteps > 1
    @assert scalerange[2] > scalerange[1] > 0

    freqs = collect(range(fmin,fmin+bandwidth,nfreqs))

    i = 1

    Scale = []; S = []

    while freqs[1] < fmax
        trace = nelderMead(booster,hist,freqs,
                    1.,1+2/booster.ndisk,0.75-1/2booster.ndisk,1-1/booster.ndisk,1e-9,
                    objective,
                    InitSimplexRegular(1e-5),
                    DefaultSimplexSampler,
                    unstuckinator;
                    maxiter=Int(1e3),
                    showtrace=false,
                    showevery=100,
                    unstuckisiter=true,)

        display(plot(freqs/1e9,getBoost1d(booster,freqs),title="new boost"))

        scale = freqs[1]/(freqs[1]+(bandwidth-overlap))
        
        freqs = collect(range(fmin+(bandwidth-overlap)*i,fmin+bandwidth*(i+1)-overlap*i,nfreqs))
        i += 1

        s = rescale(booster,hist,freqs,objective,scale,scalerange,scalesteps)

        # println("new fmin: $(freqs[1]), new fmax: $(freqs[2]), scale = $scale, s = $s")

        push!(Scale,scale); push!(S,s)
    end

    return Scale, S
end


function rescale(booster::Booster,hist::Vector{State},freqs::Array{Float64},obj::Callback,
        scale::Float64,scalerange::Tuple{Float64,Float64},scalesteps::Int)

    p_ = copy(booster.pos)
    dd = (scale-1)*pos2dist(booster.pos; disk_thickness=booster.thickness)

    b_ = 0; B = zeros(Float64,scalesteps+1)
    i_ = 0

    for i in 0:scalesteps
        move(booster,dist2pos(pos2dist(p_)+dd*Dragoon.lerp(scalerange,i/scalesteps)))
        updateHist!(booster,hist,freqs,obj)
        
        b = sum(getBoost1d(booster,freqs))
        B[i+1] = b

        if b > b_
            b_ = b
            i_ = i
        end
    end

    # display(plot(lerp.(scalerange[1],scalerange[2],(0:scalesteps)./scalesteps),B,
    #     title="rescaling, $(lerp(scalerange,i_/scalesteps))"))

    move(booster,dist2pos(pos2dist(p_)+dd*Dragoon.lerp(scalerange,i_/scalesteps)))

    return lerp(scalerange,i_/scalesteps)
end