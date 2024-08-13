
function showQuality(data,threshold=nothing)
    display(histogramB(data.obj; title="objective values",
        xlabel="objective value"))
    display(histogramB(data.opttime; title="total travel time",
        xlabel="travel time [s]"))
    display(histogramB(data.optdist; title="total travel distance",
        xlabel="trave distance [m]"))

    p1 = histogram2d(data.opttime,-data.obj; legend=false,
        xlabel="travel time [s]",ylabel="objective value")
    if !isnothing(threshold)
        hline!(p1,[-threshold])
    end
    display(p1)

    p2 = histogram2d(data.optdist,-data.obj; legend=false,
        xlabel="travel distance [m]",ylabel="objective value")
    if !isnothing(threshold)
        hline!(p2,[-threshold])
    end
    display(p2)

    p3 = histogram2d(data.optdist,data.opttime; legend=false,
        xlabel="travel distance [m]",ylabel="travel time [m]")
    display(p3)

    if !isnothing(threshold)
        sr = round(100*sum((data.obj .<= threshold))/length(data.obj); sigdigits=4)
        println("success rate: $sr")
    end

    return
end



