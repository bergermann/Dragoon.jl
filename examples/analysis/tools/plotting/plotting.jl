include("distribution.jl")
include("wiggle.jl")


function histogramB(data; kwargs...)
    p = histogram(data; legend=false,kwargs...)
    vline!(p,[minimum(data),maximum(data),sum(data)/length(data)])

    return p
end

function ridge(x,y,z=nothing; yscale=1,normalize=false,kwargs...)
    c = RGB{Float64}(0.0,0.6056031611752245,0.9786801175696073)
    a = 1

    if normalize
        m = maximum(maximum.(y))
        for i in eachindex(y)
            y[i] /= m
        end
    end

    p = plot(; size=(600,40*length(y)),legend=false,kwargs...)

    yticks!(p,1:length(y))

    if eltype(x) <: Number
        for i in reverse(eachindex(y))
            plot!(p,x,yscale*y[i].+i; fillrange=i*ones(length(y[i])),label="",c=:black,
                fillalpha=a,fillcolor=c)
            # plot!(p,x,i*ones(length(y[i])),label="";c=:black)
        end
    else
        @assert length(x) == length(y) "Arrays for x and y need same amount of entries."

        for i in reverse(eachindex(y))
            plot!(p,x[i],yscale*y[i].+i; fillrange=i*ones(length(y[i])),label="",c=:black,
                fillalpha=a,fillcolor=c)
        end
    end

    return p
end

