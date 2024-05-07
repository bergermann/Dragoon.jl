
function histogramB(data,bins=:auto)
    p = histogram(data; bins=bins,legend=false)
    vline!(p,[minimum(data),maximum(data),sum(data)/length(data)])

    return p
end