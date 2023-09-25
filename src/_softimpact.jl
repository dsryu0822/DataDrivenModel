include("../src/factorio.jl")
include("../src/DDM.jl")
include("../src/visual.jl")
const _dt = 10^(-5)

d_range = 0.1:0.0001:0.3
plan = DataFrame(idx=eachindex(d_range), d=d_range)

cd("//155.230.155.221/ty/DS");
pwd()
## -------------- begin bifurcation diagram -------------- ##
# xdots = Float64[]; ydots = Float64[]
# for dr in ProgressBar(eachrow(schedule))
#     data = factory_soft(dr.idx, dr.d)
#     idx = [false; diff(abs.(data.u) .> (dr.d/2)) .> 0]
#     sampled = data.v[idx]
#     append!(xdots, fill(dr.d, length(sampled)))
#     append!(ydots, sampled)
# end
# @time a1 = scatter(xdots, ydots,
# xlabel = L"d", ylabel = "Impact velocity",
# label = :none, msw = 0, color = :black, ms = 0.5, alpha = 0.5, size = (700, 300))
# png(a1, "soft_bifurcation")
## -------------- end bifurcation diagram -------------- ##

dr = eachrow(plan)[1]
@time data = factory_soft(dr.idx, dr.d)
look(data)
valnames = ["t" "u" "v" "cos(t)" "cos(u)" "cos(v)" "abs(t)" "abs(u)" "abs(v)" "sign(t)" "sign(u)" "sign(v)"] |> vec

_data = data[1:10:end, :]

using ProgressBars
using Random

n = nrow(_data)
sampled = rand(1:n, 12) # 230920 샘플링을 라이브러리 수와 똑같이 맞춰버리면 SingularException이 발생할 수 있음
STLSQ(_data[sampled, :], [:du, :dv], [:t, :u, :v], f_ = [cospi, abs, sign], verbose=true)


stranger = Int64[]
error_ = Float64[]
for k in ProgressBar(1:n)
    if k ∈ sampled
        continue
    end
    sampledk = [sampled; k]

    result = STLSQ(_data[sampledk, :], [:du, :dv], [:t, :u, :v], f_ = [cospi, abs, sign])
    error = result.MSE
    push!(error_, result.MSE)
    if error > eps(Float64)
        push!(stranger, k)
    end
end
scatter(log10.(error_))

subsystem = ones(Int64, n);
subsystem[stranger] .= 2;
_data[:, :subsystem] = subsystem;

scatter(_data.u, _data.v, _data.dv,
    xlabel=L"u", ylabel=L"v", zlabel=L"dv",
    color=subsystem, label=:none,
    ms=1, alpha=0.5,
    size=(800, 800))
png("temp 1")

gdf_ = groupby(_data, :subsystem)
STLSQ_ = [STLSQ(gdf, [:du, :dv], [:t, :u, :v], f_ = [cospi, abs, sign]) for gdf in gdf_]
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
Dtree = DecisionTreeClassifier(max_depth=3)
fit!(Dtree, Matrix(_data[:, [:t, :u, :v]]), subsystem)
print_tree(Dtree, 5)

x_ = [collect(data[1, [:t, :u, :v]])]
x = x_[end]

for t in ProgressBar(x[1]:_dt:50)
    s = predict(Dtree, x_[end])
    x, dx = RK4(g, s, x_[end], _dt)
    push!(x_, x)
end
_x_ = stack(x_)

uv1 = plot(data.t[1:10:end], data.u[1:10:end], label="data")
plot!(uv1, data.t[1:10:end], _x_[2, 1:10:end], color=:red, style=:dash, label="predicted")
title!(uv1, "Soft impact")

uv2 = plot(data.t[1:10:end], data.v[1:10:end], label="data")
plot!(uv2, data.t[1:10:end], _x_[3, 1:10:end], color=:red, style=:dash, label="predicted")
xlabel!(uv2, L"t")

plot(uv1, uv2, layout=(2, 1), size=(800, 800));
png("temp 2");