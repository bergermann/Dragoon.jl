using FFTW
using JLD2
using Plots



function idealBandPass(y::AbstractVector,x::AbstractVector,lld::Real=-Inf64,uld::Real=Inf64)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    y_ = copy(y)
    for i in eachindex(y)
        y_[i] *= lld < x[i] < uld
    end

    return y_
end

function idealBandPass!(y::AbstractVector,x::AbstractVector,lld::Real=-Inf64,uld::Real=Inf64)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    for i in eachindex(y)
        y[i] *= lld < x[i] < uld
    end

    return y
end

function idealBandPass(y::AbstractVector,lld::Real=0,uld::Real=1)
    @assert 0 <= lld <= 1 && 0 <= uld <= 1 "Discriminators need to be within range [0,1]"
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    y_ = copy(y)
    l = length(y)+1
    for i in eachindex(y)
        y_[i] *= lld < i/l < uld
    end

    return y_
end

function idealBandPass!(y::AbstractVector,lld::Real=0,uld::Real=1)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    l = length(y)+1
    for i in eachindex(y)
        y[i] *= lld < i/l < uld
    end

    return y
end



function idealBandStop(y::AbstractVector,x::AbstractVector,lld::Real=-Inf64,uld::Real=Inf64)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    y_ = copy(y)
    for i in eachindex(y)
        y_[i] *= !(lld < x[i] < uld)
    end

    return y_
end

function idealBandStop!(y::AbstractVector,x::AbstractVector,lld::Real=-Inf64,uld::Real=Inf64)
    @assert length(x) == length(y) "x and y require same length."
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    for i in eachindex(y)
        y[i] *= !(lld < x[i] < uld)
    end

    return y
end

function idealBandStop(y::AbstractVector,lld::Real=0,uld::Real=1)
    @assert 0 <= lld <= 1 && 0 <= uld <= 1 "Discriminators need to be within range [0,1]"
    @assert lld < uld "Lower level discriminator needs to lower than upper level discriminator."

    y_ = copy(y)
    l = i+1
    for i in eachindex(y)
        y_[i] *= !(lld < i/l < uld)
    end

    return y_
end

function idealBandStop!(y::AbstractVector,lld::Real=0,uld::Real=1)
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

function timeGate(y::AbstractVector,f::AbstractVector,llg::Real=-Inf64,ulg::Real=Inf64)
    yt, t = freq2time(y,f)
    idealBandPass!(yt,t,llg,ulg)
    yf, _ = time2freq(yt,t)

    return yf
end



@load "C:\\Users\\domin\\OneDrive\\Desktop\\testdata_3_7_sleep1.jld2"

s = R[1]
f = collect(range(18,22,length(s)))*1e9
p = plot(f/1e9,abs.(s),xlabel="f [GHz]")


yt,t = freq2time(s,f)
plot(t/1e-9,abs.(yt),xlim=[-0,50],xlabel="t [ns]")

yf,f_ = time2freq(yt,t)
plot(f_/1e9,abs.(yf),xlabel="f [GHz]")



yf = timeGate(s,f,10e-9,20e-9)
plot(p,f/1e9,abs.(yf))



# F = fftshift(fft(s))
# freqs = fftshift(fftfreq(length(x),1/(x[2]-x[1])))

# plot(x,abs.(s))
# plot!(x,real.(s))
# plot!(x,imag.(s))

# plot(freqs,abs.(F),xlim=[0,5]); vline!([1.5,2.5,3.5])

# F_ = idealBandPass(F,freqs,0.9*1.5,1.1*1.5)
# plot(freqs,abs.(F_),xlim=[0,5]); vline!([0.9*1.5,1.1*1.5])

# s_ = ifft(ifftshift(F_))

# plot(x,abs.(s_))
# plot!(x,real.(s_))
# plot!(x,imag.(s_))