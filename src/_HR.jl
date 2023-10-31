include("../src/factorio.jl")
include("../src/DDM.jl")
include("../src/visual.jl")

data = factory_HR(0, 1, (0,200))
plot(data.x, data.z, legend = :none, xlabel = "x", ylabel = "z")
hline!([1, -1], color = :red)
plot(data.t, data.z, legend = :none, xlabel = "t", ylabel = "dz")
hline!([1, -1], color = :red)

foo = vec(sum(abs2, Matrix(data[:, [:x, :y, :z]]), dims = 2))

plot(abs.((foo ./ circshift(foo, 1)) .- 1)[Not(1)], xlabel = "t", ylabel = "r", title = "algorithm 1")
bar = abs.((foo ./ circshift(foo, 1)) .- 1)

plot(data.z, legend = :none, xlabel = "t_k", ylabel = "z", title = "algorithm 1")
hline!([1, -1], color = :red)
vline!(findall(bar .> 0.1), alpha = 0.5, color = :blue, style = :dash)

cd("//155.230.155.221/ty/DS");
pwd()

look(data)
_data = data[1:end, :]; n = nrow(_data)
# plot(data.V, data.I)

using ProgressBars

normeddf = sum.(abs2, eachrow(diff(Matrix(_data[:, [:dx, :dy, :dz]]), dims = 1)))

# plot(abs.(diff(normeddf)))
# plot(normeddf)

# plot(_data.z, xlims = (0, 50000))
# hline!([1, -1], color = :blue)

jumpT = [1; findall(normeddf .> 0.02)]
vline!(jumpT, color = :red)

T_ = UnitRange.(jumpT, circshift(jumpT, -1)); pop!(T_); T_
idx_T_ = [argmax(length.(T_))]

STLSQ_ = []
for i in axes(T_, 1)
    sampled = []
    for j in 1:100
        sampled = rand(reduce(vcat, T_[[idx_T_; i]]), 200)
        X = Θ(_data[sampled, [:t, :x, :y, :z]], K = 5)
        if rank(X) ≥ max(size(X, 2))
            break
        end
    end
    if rank(X) < max(size(X, 2))
        println("sampling failed: rank = ", rank(X), " < ", max(size(X, 2)))
    end
    
    STLSQed = STLSQ(_data[sampled, :], [:dx, :dy, :dz], [:t, :x, :y, :z], K = 5, λ = 0.001)
    push!(STLSQ_, STLSQed)
    println("\ni=$i, MSE=$(STLSQed.MSE) with ", length(STLSQed.matrix.nzval))
end

STLSQ_[2]
argmin(length.(getproperty.(getproperty.(STLSQ_[Not(8)], :matrix), :nzval)))
Matrix.(getproperty.(STLSQ_[[2, 8]], :matrix))[1]
Matrix.(getproperty.(STLSQ_[[4, 8]], :matrix))[2]

Θ(_data[1:10, [:t, :x, :y, :z]], K = 5)

sampled = rand(reduce(vcat, T_[[idx_T_; 2]]), 200)
STLSQed = STLSQ(_data[sampled, :], [:dx, :dy, :dz], [:t, :x, :y, :z], K = 5, λ = 0.001)
α = Matrix(STLSQ_[4].matrix)

stranger = Int64[]
error_ = Float64[]
for k in ProgressBar(1:n)
    error = sum(abs2, collect(_data[k, [:dx, :dy, :dz]])' .- Θ(_data[k, [:t, :x, :y, :z]], K = 3, f_ = [cos])*α)
    push!(error_, error)
    if error > 1e-8
        push!(stranger, k)
    end
end
scatter(log10.(error_))

subsystem = ones(Int64, n);
subsystem[stranger] .= 2;
_data[:, :subsystem] = subsystem;


idx_left = findall([count(stranger .∈ Ref(TT)) for TT in T_] .> 100)
T_ = T_[idx_left]

idx_T_ = [argmax(length.(T_))]

STLSQ_ = []
for i in axes(T_, 1)
    sampled = []
    for j in 1:100
        sampled = rand(reduce(vcat, T_[[idx_T_; i]]), 50)
        if rank(Θ(_data[sampled, [:t, :x, :y, :z]], K = 3, f_ = [cos])) ≥ 39
            break
        end
    end
    
    STLSQed = STLSQ(_data[sampled, :], [:dx, :dy, :dz], [:t, :x, :y, :z], K = 3, f_ = [cos], λ = 0.001)
    push!(STLSQ_, STLSQed)
    println("\ni=$i, MSE=$(STLSQed.MSE) with ", length(STLSQed.matrix.nzval))
end

α = Matrix(STLSQ_[3].matrix)

stranger = Int64[]
error_ = Float64[]
for k in ProgressBar(1:n)
    error = sum(abs2, collect(_data[k, [:dx, :dy, :dz]])' .- Θ(_data[k, [:t, :x, :y, :z]], K = 3, f_ = [cos])*α)
    push!(error_, error)
    if error > 1e-8
        push!(stranger, k)
    end
end
scatter(log10.(error_))
_data.subsystem[findall(subsystem .> 1) ∩ stranger] .= 3;


plot(_data.x[1:10:end], _data.z[1:10:end],
    xlabel=L"x", ylabel=L"z",
    color=_data.subsystem[1:10:end], label=:none,
    size=(800, 800))
hline!([1, -1], color = :black, style = :dash)
png("temp 1")

gdf_ = groupby(_data, :subsystem)
STLSQ_ = [STLSQ(gdf, [:dx, :dy, :dz], [:t, :x, :y, :z], K = 3, f_ = [cos]) for gdf in gdf_]
Matrix.(getproperty.(STLSQ_, :matrix))[1]
Matrix.(getproperty.(STLSQ_, :matrix))[2]
Matrix.(getproperty.(STLSQ_, :matrix))[3]

function factory_STLSQ(STLSQed)
    function f(s, x)
        return [1; vec(Θ(x, K = 3, f_ = [cos]) * STLSQed[s].matrix)]
    end
    return f
end
g = factory_STLSQ(STLSQ_)


using DecisionTree
my_depth = 3
Dtree = DecisionTreeClassifier(max_depth=my_depth)
features = Matrix(_data[:, [:t, :x, :y, :z]])
fit!(Dtree, features, _data.subsystem)
print_tree(Dtree, my_depth)
acc = sum(_data.subsystem .== predict(Dtree, features)) / length(subsystem)
println(acc)

x_ = [collect(data[1, [:t, :x, :y, :z]])]
x = x_[end]

dt = 10^(-3)
for t in ProgressBar(0:dt:200)
    s = predict(Dtree, x_[end])
    x, dx = RK4(g, s, x_[end], dt)
    push!(x_, x)
end
_x_ = stack(x_)[:, 1:(end-2)]

uv1 = plot(data.t[1:10:end], data.x[1:10:end], ylabel = L"x", label="data")
plot!(uv1, data.t[1:10:end], _x_[2, 1:10:end], color=:red, alpha = 0.5, label="predicted")
title!(uv1, "Hindmarsh-Rose neuronal model")

uv2 = plot(data.t[1:10:end], data.z[1:10:end], ylabel = L"z", label="data")
plot!(uv2, data.t[1:10:end], _x_[4, 1:10:end], color=:red, alpha= 0.5, label="predicted")
xlabel!(uv2, L"t")
plot(uv1, uv2, layout=(2, 1), size=(800, 800));
png("temp 2");

# ---

@time data = factory_buck(dr.idx, dr.E, (0.495, 0.55))
x_ = [collect(data[1, [:V, :I]])]
x = x_[end]

for ramp in ProgressBar(data.Vr)
    s = predict(Dtree, [x_[end]; ramp])
    x, dx = RK4(g, s, x_[end], dt)
    push!(x_, x)
end
_x_ = stack(x_)

uv1 = plot(data.t[1:10:end], data.V[1:10:end], ylabel = L"V", label="data")
plot!(uv1, data.t[1:10:end], _x_[1, 1:10:end], color=:red, style=:dash, label="predicted")
title!(uv1, "Buck converter")

uv2 = plot(data.t[1:10:end], data.I[1:10:end], ylabel = L"I", label="data")
plot!(uv2, data.t[1:10:end], _x_[2, 1:10:end], color=:red, style=:dash, label="predicted")
xlabel!(uv2, L"t")
plot(uv1, uv2, layout=(2, 1), size=(1600, 900))