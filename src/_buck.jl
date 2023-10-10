include("../src/factorio.jl")
include("../src/DDM.jl")
include("../src/visual.jl")
const dt = 10^(-7)

E_range = 15:0.001:40
plan = DataFrame(idx=eachindex(E_range), E=E_range)

cd("//155.230.155.221/ty/DS");
pwd()

dr = eachrow(plan)[end]
@time data = factory_buck(dr.idx, dr.E)
look(data)
_data = data[1:end, :]
# plot(data.V, data.I)

using ProgressBars
using Random

n = nrow(_data)
sampled = rand(1:n, 3) # 230920 샘플링을 라이브러리 수와 똑같이 맞춰버리면 SingularException이 발생할 수 있음
α = STLSQ(_data[sampled, :], [:dV, :dI], [:V, :I], verbose=true).matrix

stranger = Int64[]
error_ = Float64[]
for k in ProgressBar(1:n)
    error = sum(abs2, collect(_data[k, [:dV, :dI]])' .- [0; collect(_data[k, [:V, :I]])]'α)
    push!(error_, error)
    if error > 1e-8
        push!(stranger, k)
    end
end
scatter(log10.(error_))

subsystem = ones(Int64, n);
subsystem[stranger] .= 2;
_data[:, :subsystem] = subsystem;

plot(_data.V, _data.I,
    xlabel=L"V", ylabel=L"I",
    color=subsystem, label=:none,
    ms=1, alpha=0.5,
    size=(800, 800))
png("temp 1")

gdf_ = groupby(_data, :subsystem)
STLSQ_ = [STLSQ(gdf, [:dV, :dI], [:V, :I]) for gdf in gdf_]
STLSQ_[1]
STLSQ_[2]

function factory_STLSQ(STLSQed)
    function f(s, x)
        return vec(poly_basis(x, 1)' * STLSQed[s].matrix)
    end
    return f
end
g = factory_STLSQ(STLSQ_)


using DecisionTree
my_depth = 8
Dtree = DecisionTreeClassifier(max_depth=my_depth, n_subfeatures=2)
features = Matrix(_data[:, [:V, :I, :Vr]])
fit!(Dtree, features, subsystem)
print_tree(Dtree, my_depth)
acc = sum(subsystem .== predict(Dtree, features)) / length(subsystem)
println(acc)

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