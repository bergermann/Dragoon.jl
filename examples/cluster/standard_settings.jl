
mutable struct Settings
    f0::AbstractFloat
    df::AbstractFloat
    nf::Int
    ndisk::Int
    eps::AbstractFloat
    tand::AbstractFloat
end


s = Settings(
    22.025e9,   # f0::AbstractFloat
    50e6,       # df::AbstractFloat
    10,         # nf::Int
    20,         # ndisk::Int
    24.0,       # eps::AbstractFloat
    0.0,        # tand::AbstractFloat
)