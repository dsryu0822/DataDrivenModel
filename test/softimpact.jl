include("../src/factorio.jl")
include("../src/ML.jl")
include("../src/DDM.jl")


using Plots, LaTeXStrings; default(msw = 0, color = :black)

d_range = 0.1:0.0001:0.3
plan = DataFrame(idx = eachindex(d_range), d = d_range)

cd("//155.230.155.221/ty/DS")

dr = eachrow(plan)[1]
data = factory_soft(dr.idx, dr.d)
valnames = ["t" "u" "v" "cos(t)" "cos(u)" "cos(v)" "abs(t)" "abs(u)" "abs(v)" "sign(t)" "sign(u)" "sign(v)"] |> vec

_data = data[1:10:end, :]
# points = Matrix(_data[:, [:u, :v, :dv]] |> col_normalize)'


# kmeansed = kmeans(log10.(abs.(_data.dv))', 2); subsystem = kmeansed.assignments
# @time dbscaned = dbscan(points, 0.0025); subsystem = dbscaned.assignments;
# subsystem = 1 .+ (abs.(_data.dv) .> 1)
# _data[:, :subsystem]  = subsystem;

using ProgressBars
using Random

n = nrow(_data)
sampled = rand(1:n, 12) # 230920 샘플링을 라이브러리 수와 똑같이 맞춰버리면 SingularException이 발생할 수 있음
X = col_func(Matrix(_data[sampled, [:t, :u, :v]]), [cospi, abs, sign])
Y = Matrix(_data[sampled, [:du, :dv]])
STLSQ(X, Y, λ = 0.1)

stranger = Int64[]
error_ = Float64[]
for k in ProgressBar(1:n)
    if k ∈ sampled continue end
    sampledk = [sampled; k]
    X = col_func(Matrix(_data[sampledk, [:t, :u, :v]]), [cospi, abs, sign])
    Y = Matrix(_data[sampledk, [:du, :dv]])
    try
        result = STLSQ(X, Y, λ = 0.1, verbose = false)
        error = result.MSE
        push!(error_, result.MSE)
        if log10(error) > -20
            push!(stranger, k)
        end
    catch err
        if err isa SingularException
            push!(error_, 10^(-30))
        end
    end
end
scatter(log10.(error_))

# Z = col_func(Matrix(_data[[sampled; 1], [:t, :u, :v]]), [cospi, abs, sign])
# det(Z)
# Z[:,1] .== Z[:,7]


subsystem = ones(Int64, n); subsystem[stranger] .= 2;
_data[:, :subsystem]  = subsystem;

scatter(_data.u, _data.v, _data.dv,
xlabel = L"u", ylabel = L"v", zlabel = L"dv",
color = subsystem, label = :none,
ms = 1, alpha = 0.5,
size = (800, 800))
png("230920_1")

gdf_ = groupby(_data, :subsystem)
STLSQ_ = []
for gdf in gdf_
    # X = poly_basis(Matrix(gdf[:, [:t, :u, :v]]), 3)
    X = col_func(Matrix(gdf[:, [:t, :u, :v]]), [cospi, abs, sign])
    Y = Matrix(gdf[:, [:du, :dv]])
    push!(STLSQ_, STLSQ(X, Y, λ = 0.1))
end
STLSQ_[1]
STLSQ_[2]

function factory_STLSQ(STLSQed)
    function f(s, x)
        return [1; vec(col_func(x', [cospi, abs, sign]) * STLSQed[s].matrix)]
    end
    return f
end
g = factory_STLSQ(STLSQ_)


using DecisionTree
Dtree = DecisionTreeClassifier(max_depth = 3)
fit!(Dtree, Matrix(_data[:, [:t, :u, :v]]), subsystem)
print_tree(Dtree, 5)

x_ = [collect(data[1,1:3])]
x = x_[end]

dt = 10^(-5)
for t in ProgressBar(x[1]:dt:50)
    s = predict(Dtree, x_[end])
    x, dx = RK4(g, s, x_[end], dt)
    push!(x_, x)
    push!(dx_, dx)
end

uv1 = plot(data.t[1:10:end], data.u[1:10:end], label = "data")
plot!(uv1, data.t[1:10:end], stack(x_)[2,1:10:end], color = :red, style = :dash, label = "predicted")

uv2 = plot(data.t[1:10:end], data.v[1:10:end], label = "data")
plot!(uv2, data.t[1:10:end], stack(x_)[3,1:10:end], color = :red, style = :dash, label = "predicted")

plot(uv1, uv2, layout = (2,1), size = (800, 800)); png("230920_2")