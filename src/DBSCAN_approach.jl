@time using Plots; mm = Plots.mm
using NearestNeighbors
include("nonsmooth.jl")

function patrol(v)
    x, s = v

    if s > 0
        return [-0.1exp(abs(x)), 0]
    elseif s < 0
        return [0.1exp(abs(x)), 0]
    end
end

v_ = [Float64[3,1]]
d_ = [patrol(v_[1])]
dt = 0.01
for t in 0:dt:100
    push!(v_, RK4(patrol, v_[end], dt))
    if v_[end][1] ≥ 1
        v_[end][2] = 1
    elseif v_[end][1] ≤ -1
        v_[end][2] = -1
    end
    push!(d_, patrol(v_[end]))
end
X = stack(v_)'[:,1]
Y = stack(d_)'[:,1]
datargs = (; xformatter = (x -> x*dt), xlabel = "t", label = "Data",
            color = :black, legend = :none)
p1 = plot(X[:,1]; datargs..., title = "Time Evolution", ylabel = "x")
p2 = plot(Y[:,1]; datargs..., ylabel = "ẋ")
p3 = scatter(X,Y, xlabel = "x", ylabel = "dx/dt", label = "(x, ẋ)", color = :blacK,  title = "Scatter Plot: x vs ẋ")
plot(p1, p2, p3; margin = 5mm, size = (900, 600), layout = @layout [[a;b] c])

using Clustering

XY = [X Y]'
rslt_dbscn = dbscan(XY, 0.1)
# rslt_dbscn |> propertynames
rslt_dbscn.clusters[2].core_indices
plot(
scatter(X,Y, xlabel = "x", ylabel = "dx/dt", label = "(x, ẋ)", color = rslt_dbscn.assignments, msw = 0,  title = "Scatter Plot: x vs ẋ", legend = :none),
plot(X[:,1]; datargs..., color = rslt_dbscn.assignments, title = "Time Evolution", lw = 2),
)
plot(X[:,1]; datargs..., color = rslt_dbscn.assignments, title = "Time Evolution", lw = 2, size = (1000, 300))
plot(Y[:,1]; datargs..., color = rslt_dbscn.assignments, title = "Time Evolution", lw = 2, ylabel = "ẋ")
png("50.png")
θ₁, θ₂ = extrema(X[argjump(Y)])


bit_subsystem = spzeros(Bool, length(X), 2)
bit_subsystem[rslt_dbscn.clusters[1].core_indices, 1] .= 1
bit_subsystem[rslt_dbscn.clusters[2].core_indices, 2] .= 1

Ξ = STLSQ(col_prod(bit_subsystem, poly_basis(reshape(X, :, 1),4)), Y, 0.1)

# KDTreeXY = KDTree(XY)
function foo(v)
    x,s = v
    subsystem = ifelse(s > 0, [1, 0], [0, 1])
    z = col_prod(subsystem, poly_basis([x], 4))
    return [(z * Ξ)..., 0]
end

v_ = [Float64[3,1]]
d_ = [Float64[-2,0]]
dt = 0.01
for t in 0:dt:100
    push!(v_, RK4(foo, [v_[end][1]; v_[end][end]], dt))
    if v_[end][1] ≥ θ₂
        v_[end][end] = 1
    elseif v_[end][1] ≤ θ₁
        v_[end][end] = -1
    end
    push!(d_, foo(v_[end]))
end
myX = stack(v_)'[:,1]
myY = stack(d_)'[:,1]
plot(X; datargs..., legend = :topright)
prdargs = (; xformatter = (x -> x*dt), xlabel = "t", label = "Predicted",
            color = :blue, ls = :dash)
plot!(myX; prdargs...)