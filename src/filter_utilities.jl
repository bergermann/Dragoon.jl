using FFTW


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
    t = fftshift(fftfreq(length(f),1/(f[2]-f[1])))

    return yt, t
end

function time2freq(yt::AbstractArray,t::AbstractArray)
    yf = fftshift(fft(yt))
    f = ifftshift(fftfreq(length(t),1/(t[2]-t[1])))

    return yf, f
end

function timeGate(y::AbstractVector,f::AbstractVector,llg::Real=-Inf64,ulg::Real=Inf64)
    yt, t = freq2time(y,f)
    idealBandPass!(yt,t,llg,ulg)
    yf = time2freq(yt,t)

    return y
end



f = collect(18:0.001:22)*1e9
a = 1.0*cispi.(2*0.03e-9*x)
b = 0.1*cispi.(2*2e-9*x)
c = 0.05*cispi.(2*3e-9*x)
d = 0.01*(2*rand(length(x)).-1)

s = a+b+c+d

plot(f/1e9,real.(s))

yt,t = freq2time(s,f)

plot(t,abs.(yt),xlim=[-10,10]/1e9)


yf,f_ = time2freq()









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