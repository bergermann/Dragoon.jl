using FFTW
using JLD2
using Plots



function bandPass(y::AbstractVector,x::AbstractVector,lld::Real=-Inf64,uld::Real=Inf64; ramp::Real=0)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."
    @assert ramp <= uld-lld "Ramp needs to be smaller than uld-lld."

    y_ = copy(y)
    if ramp == 0
        @simd for i in eachindex(y)
            y_[i] *= lld < x[i] < uld
        end
    else
        @simd for i in eachindex(y)
            y_[i] *= clamp((x[i]-lld)/2ramp+0.5,0.,1.)*clamp(-(x[i]-uld)/2ramp+0.5,0.,1.)
        end
    end

    return y_
end

function bandPass!(y::AbstractVector,x::AbstractVector,lld::Real=-Inf64,uld::Real=Inf64; ramp::Real=0)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."
    @assert ramp <= uld-lld "Ramp needs to be smaller than uld-lld."

    if ramp == 0
        @simd for i in eachindex(y)
            y[i] *= lld < x[i] < uld
        end
    else
        @simd for i in eachindex(y)
            y[i] *= clamp((x[i]-lld)/2ramp+0.5,0.,1.)*clamp(-(x[i]-uld)/2ramp+0.5,0.,1.)
        end
    end

    return y
end

function bandPass(y::AbstractVector,lld::Real=0,uld::Real=1; ramp::Real=0)
    @assert 0 <= lld <= 1 && 0 <= uld <= 1 "Discriminators need to be within range [0,1]"
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."
    @assert ramp <= uld-lld "Ramp needs to be smaller than uld-lld."

    y_ = copy(y)
    l = length(y)+1
    if ramp == 0
        @simd for i in eachindex(y)
            y_[i] *= lld < i/l < uld
        end
    else
        @simd for i in eachindex(y)
            y_[i] *= clamp((i/l-lld)/2ramp+0.5,0.,1.)*clamp(-(i/l-uld)/2ramp+0.5,0.,1.)
        end
    end

    return y_
end

function bandPass!(y::AbstractVector,lld::Real=0,uld::Real=1; ramp::Real=0)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."
    @assert ramp <= uld-lld "Ramp needs to be smaller than uld-lld."

    l = length(y)+1
    if ramp == 0
        @simd for i in eachindex(y)
            y[i] *= lld < i/l < uld
        end
    else
        @simd for i in eachindex(y)
            y[i] *= clamp((i/l-lld)/2ramp+0.5,0.,1.)*clamp(-(i/l-uld)/2ramp+0.5,0.,1.)
        end
    end

    return y
end



function bandStopIdeal(y::AbstractVector,x::AbstractVector,lld::Real=-Inf64,uld::Real=Inf64)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    y_ = copy(y)
    for i in eachindex(y)
        y_[i] *= !(lld < x[i] < uld)
    end

    return y_
end

function bandStopIdeal!(y::AbstractVector,x::AbstractVector,lld::Real=-Inf64,uld::Real=Inf64)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    for i in eachindex(y)
        y[i] *= !(lld < x[i] < uld)
    end

    return y
end

function bandStopIdeal(y::AbstractVector,lld::Real=0,uld::Real=1)
    @assert 0 <= lld <= 1 && 0 <= uld <= 1 "Discriminators need to be within range [0,1]"
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    y_ = copy(y)
    l = i+1
    for i in eachindex(y)
        y_[i] *= !(lld < i/l < uld)
    end

    return y_
end

function bandStopIdeal!(y::AbstractVector,lld::Real=0,uld::Real=1)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    l = i+1
    for i in eachindex(y)
        y[i] *= !(lld < i/l < uld)
    end

    return y
end

function freq2time(yf::AbstractArray,f::AbstractArray)
    yt = ifftshift(ifft(yf))
    t = fftshift(fftfreq(length(f),(length(f)-1)/(f[end]-f[1])))

    return yt, t
end

function time2freq(yt::AbstractArray,t::AbstractArray)
    yf = fft(yt)
    dt = (length(t)-1)/(t[end]-t[1])
    f = fftshift(fftfreq(length(t),dt)).+dt

    return yf, f
end

function timeGate(y::AbstractVector,f::AbstractVector,llg::Real=-Inf64,ulg::Real=Inf64; ramp::Real=0.)
    yt, t = freq2time(y,f)
    bandPass!(yt,t,llg,ulg; ramp=ramp)
    yf, _ = time2freq(yt,t)

    return yf
end


# @load "\\\\inst3\\data\\Benutzer\\bergermann\\Desktop\\testdata_3_7_sleep1.jld2"
# @load "C:\\Users\\domin\\OneDrive\\Desktop\\testdata_3_7_sleep1.jld2"

# s = R[1]
# f = collect(range(18,22,length(s)))*1e9
# p = plot(f/1e9,abs.(s),xlabel="f [GHz]")


# yt,t = freq2time(s,f)
# plot(t/1e-9,abs.(yt),xlim=[-20,50],xlabel="t [ns]")

# yf,f_ = time2freq(yt,t)
# plot(f_/1e9,abs.(yf),xlabel="f [GHz]")



# yf = timeGate(s,f,10e-9,20e-9; ramp=3e-9)
# plot(p,f/1e9,abs.(yf))



