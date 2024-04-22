
mutable struct Settings
    f0::Float64
    df::Float64
    nf::Int
    ndisk::Int
    eps::Float64
    tand::Float64
end


s = Settings(
    22.025e9,   # f0::Float64
    50e6,       # df::Float64
    10,         # nf::Int
    20,         # ndisk::Int
    24.0,       # eps::Float64
    0.0,        # tand::Float64
)


import Base: println
function Base.println(s::Settings)
    println("--- Settings ---")
    println("Center frequency: $(s.f0)")
    println("Span frequency:   $(s.df)")
    println("Frequency points: $(s.nf)")
    println("Disc amount:      $(s.ndisk)")
    println("Disc epsilon:     $(s.eps)")
    println("Disc loss (tand): $(s.tand)")
end