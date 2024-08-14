
function showQuality(data,threshold=nothing)
    display(histogramB(data.obj; title="Objective Values",
        xlabel="Objective Value"))
    display(histogramB(data.opttime; title="Total Travel Time",
        xlabel="Travel Time [s]"))
    display(histogramB(data.optdist; title="Total Travel Distance",
        xlabel="Travel Distance [m]"))

    p1 = histogram2d(data.opttime,-data.obj; label="",
        xlabel="Travel Time [s]",ylabel="Objective Value")
    if !isnothing(threshold)
        hline!(p1,[-threshold])
    end
    display(p1)

    p2 = histogram2d(data.optdist,-data.obj; label="",
        xlabel="Travel Distance [m]",ylabel="Objective Value")
    if !isnothing(threshold)
        hline!(p2,[-threshold])
    end
    display(p2)

    p3 = histogram2d(data.optdist,data.opttime; legend=false,
        xlabel="Travel Distance [m]",ylabel="Travel Time [m]")
    display(p3)

    if !isnothing(threshold)
        sr = round(100*sum((data.obj .<= threshold))/length(data.obj); sigdigits=4)
        println("success rate: $sr")
    end

    return
end



