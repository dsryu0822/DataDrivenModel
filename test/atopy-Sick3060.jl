@time using CSV, DataFrames, Plots
using StatsBase

include("../src/Utils.jl")
include("../src/nonsmooth.jl")

# Oₘ = CSV.read("G:/DDM/atopy.txt", DataFrame, delim = '\t', header = ["P","B","D","R","K","G","kappa","t"], skipto = 10000)
Oₘ = CSV.read("G:/DDM/atopy.txt", DataFrame, delim = '\t', header = ["P","B","D","R","K","G","kappa","t"])
dt = Oₘ.t[2] - Oₘ.t[1]
# sP = sqrt(sum(abs2, Oₘ.P))
# sB = sqrt(sum(abs2, Oₘ.B))
# sD = sqrt(sum(abs2, Oₘ.D))
# Oₘ.P /= sP
# Oₘ.B /= sB
# Oₘ.D /= sD

X = Matrix(Oₘ[1:1:(10000-1),1:3])
Y = ((Matrix(Oₘ[2:1:10000,1:3]) - Matrix(Oₘ[1:1:(10000-1),1:3])) / dt)

a1 = plot(X[:,3], color = :black, label = :none)
scatter!(a1, argjump(Y[:,3]), X[argjump(Y[:,3]),3], label = "Jump", color = :white)

plot(X[:,1])
plot(X[:,2])

using Clustering

# dbs |> propertynames
dbs = dbscan([X Y]', 12)
dbs.clusters
plot(
    plot(X[:,1], ylabel = "P", lw = 2, legend = :none, color = dbs.assignments),
    plot(X[:,2], ylabel = "B", lw = 2, legend = :none, color = dbs.assignments),
    plot(X[:,3], ylabel = "D", lw = 2, legend = :none, color = dbs.assignments),
    a1,
    xlabel = "time step", size = 80 .* (12, 9), margin = 5Plots.mm
) # png("61.png")

θ₁, θ₂ = extrema(X[argjump(Y[:,1]),1])
nc = length(dbs.clusters)

bit_subsystem = spzeros(Bool, size(Y,1), nc)
bit_subsystem[dbs.clusters[1].core_indices, 1] .= 1
bit_subsystem[dbs.clusters[2].core_indices, 2] .= 1
bit_subsystem[dbs.clusters[3].core_indices, 3] .= 1
bit_subsystem[dbs.clusters[4].core_indices, 4] .= 1
Inc = I(nc)

Θ = col_prod(bit_subsystem, poly_basis(X,4))
Ξ = STLSQ(Θ, Y, 0)
function foo(v)
    s = Int(v[end])
    subsystem = Inc[s,:]
    z = col_prod(subsystem, poly_basis(v[1:(end-1)], 4))
    return [(z * Ξ)..., 0]
end

v_ = [[X[1,:]..., 1]]
for t in 2:9999
    push!(v_, RK4(foo, v_[end], dt))
    if v_[end][1] ≥ θ₂
        v_[end][end] = 4
    elseif v_[end][1] ≤ θ₁
        v_[end][end] = 3
    end
end
# for t in 1:10000
#     push!(v_, RK4(foo, v_[end], dt))
#     if v_[end][1] ≥ θ₂
#         v_[end][end] = 2
#     elseif v_[end][1] ≤ θ₁
#         v_[end][end] = 1
#     end
# end
myX = stack(v_)'
plot(X[:,3], color = :black, ylabel = "B", xlabel = "time step", legend = :none)
plot!(myX[:,3], color = :blue, ls = :dash)

[X[:,2] myX[:,2]]
plot(abs.(X[2:end,2] - myX[2:end,2]), color = :black, ylabel = "|B - B̂|", xlabel = "time step", legend = :none, yscale = :log10, title = "Absolute Error")