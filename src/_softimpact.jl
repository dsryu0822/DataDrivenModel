include("../src/factorio.jl")
include("../src/DDM.jl")
include("../src/visual.jl")
const _dt = 10^(-5)

d_range = 0.1:0.0001:0.3
plan = DataFrame(idx=eachindex(d_range), d=d_range)

dr = first(eachrow(plan))
@time data = factory_soft(DataFrame, dr.idx, dr.d, tspan = [0, 100])
look(data)
half = nrow(data) ÷ 2
cut = 10_000_000
# cut = 1_000_000
_data = data[1:cut, :]
_data = deepcopy(data)

n = nrow(_data)
sampled = rand(1:n, 15) # 230920 샘플링을 라이브러리 수와 똑같이 맞춰버리면 SingularException이 발생할 수 있음
# sampled = 1:12 # 230920 샘플링을 라이브러리 수와 똑같이 맞춰버리면 SingularException이 발생할 수 있음
f = SINDy(_data[sampled, :], [:dt, :du, :dv], [:t, :u, :v], f_ = [sign, abs, cospi])
print(f, ["t", "u", "v"])

error_ = norm.(eachrow(Matrix(_data[:, [:dt, :du, :dv]])) .- f.(eachrow(Matrix(_data[:, [:t, :u, :v]]))))
bit_alien = error_ .> 1e-5
# scatter(log10.(error_)[1:100:end])

subsystem = ones(Int64, n);
subsystem[bit_alien] .= 2;
_data[:, :subsystem] = subsystem;
plot(_data.u[1:100:end], _data.v[1:100:end], xlabel=L"u", ylabel=L"v", color=subsystem[1:100:end], label=:none, ms=1, alpha=0.5, size=(800, 800))

gdf_ = groupby(_data[1:1000:end, :], :subsystem)
f_ = [SINDy(gdf, [:dt, :du, :dv], [:t, :u, :v], f_ = [sign, cospi], λ = 1e-2) for gdf in gdf_]
print(f_[1], ["t", "u", "v"])
print(f_[2], ["t", "u", "v"])

bit_chgd = (subsystem .!= circshift(subsystem, 1)) .|| (subsystem .!= circshift(subsystem, -1)) 
__data = _data[findall(bit_chgd) ∪ rand(findall(.!bit_chgd), count(bit_chgd)), :]
# __data = deepcopy(_data)

using DecisionTree
Dtree = DecisionTreeClassifier(max_depth=3, n_subfeatures=3)
features = Matrix(__data[:, [:t, :u, :v]])
fit!(Dtree, features, __data.subsystem); print_tree(Dtree, 3)
acc = sum(__data.subsystem .== predict(Dtree, features)) / length(__data.subsystem)

dt = 1e-5
x = collect(data[end, [:t, :u, :v]])
y = factory_soft(DataFrame, dr.idx, dr.d, ic = x)
ŷ = DataFrame(solve(f_, x, dt, y.t, Dtree), ["t", "u", "v"]);
plot(xlabel = L"t", ylabel = L"u(t)", size = (800, 300), margin = 3mm)
plot!(y.t[1:100:end], y.u[1:100:end], label = "true", lw = 2)
plot!(y.t[1:100:end], ŷ.u[1:100:end], label = "pred", lw = 2, ls = :dash, color = :red)
hline!([-.05, .05], label = "bound", color = :blue)
png("112")
plot(log10.(abs.(y.u .- ŷ.u)))