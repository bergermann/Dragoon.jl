

using BoostFractor
using Plots

include("transformer_optim_utilities_old.jl")
include("transformer_utils.jl")


# ich habe die objective function so abgeändert, dass sie immer auf dem aktuellen boosterstate
# s agiert. die nötigen tilts stelle ich vorher außerhalb mit der tilts! funktion ein
# 
# ich bevorzuge wenn variablen die lokal in der funktion verwendet werden, als argumente
# an die funktion gegeben werden, daher habe ich freqs und R0 als argumente hinzugefügt
function objective_fun(s,freqs,R0)
    r = ref(freqs,s); #r_sum = sum(r)

    # da r, R0 vektoren von vektoren sind, braucht man x -> abs.(x) um abs.
    # auf jeden untervektor anzuwenden 
    return sum(sum(x->abs.(x),r - R0))
end


n_disk = 1
deg_max = 0.5

Mmax = 3; Lmax = 2

freq_min = 21.5e9; freq_max = 22.5e9; n_freq = 20
freqs = collect(range(freq_min,freq_max,n_freq));
freqsplot = collect(range(21.5e9,22.5e9,100));

d = [6.5]*1e-3

s = setup(d,Mmax,Lmax);

Bp, Rp, b, b_sump, rp = boost(freqsplot,s)
plot(freqsplot/1e9,b_sump)



tilts = collect(range(-1,1,11))*deg2rad(deg_max)
Obj = zeros(length(tilts),length(tilts))

tilt!(s.sbdry,[0],[0])
B, R, b, b_sum0, r = boost(freqs,s); R0 = r #sum(r)

objective_fun(s,freqs,R0)

for j in eachindex(tilts), i in eachindex(tilts)
    tilt!(s.sbdry,[tilts[i]],[tilts[j]])
    
    Obj[i,j] = objective_fun(s,freqs,R0)
end

contourf(tilts,tilts,Obj)





# function NelderMead(f,s,sum_deg, degs,boosts, offset = 0.0005 ; α = 1, β = 2, γ = 1/2, δ = 1/2, MaxIter = 2000, Tol = 10e-6)
#     tilt_x = s.sbdry.relative_tilt_x
#     tilt_y = s.sbdry.relative_tilt_y
#     n = length(tilt_x)
#     n_disk = length(d)
#     x0 = vcat(tilt_x,tilt_y)
#     simplex = [x0]
    
#     TiltChange = [i for i in 1:2*n_disk+1 if i != n_disk+1]

#     for i in TiltChange
#         e = zeros(2*n)
#         e[2*i] = e[1+2*i] = offset
#         push!(simplex, x0 + e)
#     end
    
#     f_dict = Dict{Vector{Float64}, Float64}()
#     n = length(simplex) - 1
#     for i in 1:MaxIter
#         println("Iteration $i")
#         # 1)
#         simplex = sort_f(simplex, f, f_dict)
        
#         # 2)
#         c = sum(simplex[1:n]) / n 

#         xr = c + α * (c - simplex[end])
#         fxr = f(xr)
#         if f(simplex[1]) <= fxr && fxr <= f(simplex[n])
#             simplex[end] = xr
        
#         # 3)
#         elseif fxr <= f(simplex[1])

#             xs = c + β * (xr - c)
#             f(xs) < fxr ? simplex[end] = xs : simplex[end] = xr 

#         # 4)
#         elseif f(simplex[n]) <= fxr && fxr <= f(simplex[end])
#             xoc = c + γ * (xr - c)
#             if f(xoc) <= fxr
#                 simplex[end] = xoc
            
#             else
#                 # 6)
#                 for j in 2:n+1
#                     simplex[j] = simplex[1] + δ * (simplex[j] - simplex[1])
#                 end
#             end
        
#         # 5)
#         else #f(xr) >= f(simplex[end])
#             xic = c - γ * (xr - c)
#             if f(xic) < fxr
#                 simplex[end] = xic
#             end

#         end

#         if i % 10 == 0
#             println(f(simplex[1]))
#             n_tarray = length(tilt_x)
#             s.sbdry.relative_tilt_x = best_tilt[1:n_tarray]; s.sbdry.relative_tilt_y = best_tilt[n_tarray+1:end]
#             B, R, b, b_sum_optim, r = boost(freqs, s)

#             plot!(p1, freqs/1e9, b_sum_optim, label = "Iteration $i")

#         end
#         push!(sum_deg, rad2deg(sum(abs.(simplex[1]))/2))
#         push!(boosts, f(simplex[1]))
#     end
#     return simplex[1], f(simplex[1])
# end

# function sort_f(v::Vector{Vector{Float64}}, f::Function, f_dict::Dict)
#     for x in v
#         haskey(f_dict, x) || (f_dict[x] = f(x))
#     end
#     ziptilt = [(x, f_dict[x]) for x in v]
#     sort!(ziptilt, by = x -> x[2])
#     simplex = [x[1] for x in ziptilt]
#     return simplex
# end

