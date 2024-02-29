include("../src/DDM.jl")
include("../src/factorio.jl")
include("../src/visual.jl")

E_range = 15:0.001:40
plan = DataFrame(idx=eachindex(E_range), E=E_range)

# cd("//155.230.155.221/ty/DS");
pwd()
dr = eachrow(plan)[end]
@time data = factory_buck(DataFrame, dr.idx, dr.E)
look(data)
half = nrow(data) ÷ 2
_data = data[1:half, :]

plot(_data.dI, xformatter = x -> x/exp10(7), legend = :none, xlabel = L"t", ylabel = L"dI")
n = nrow(_data)
# sampled = rand(1:n, 15) # 230920 샘플링을 라이브러리 수와 똑같이 맞춰버리면 SingularException이 발생할 수 있음
sampled = 1:12 # 230920 샘플링을 라이브러리 수와 똑같이 맞춰버리면 SingularException이 발생할 수 있음
f = SINDy(_data[sampled, :], [:dV, :dI], [:V, :I], verbose=true)
print(f, ["V", "I"])

error_ = norm.(eachrow(Matrix(_data[:, [:dV, :dI]])) .- f.(eachrow(Matrix(_data[:, [:V, :I]]))))
bit_alien = error_ .> 1e-4
scatter(log10.(error_)[1:100:end], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")

subsystem = ones(Int64, n);
subsystem[bit_alien] .= 2;
_data[:, :subsystem] = subsystem;
plot(_data.V[1:100:end], _data.I[1:100:end], xlabel=L"V", ylabel=L"I", color=subsystem[1:100:end], label=:none, ms=1, alpha=0.5, size=(800, 800))
__data = _data[(subsystem .!= circshift(subsystem, 1)) .|| (subsystem .!= circshift(subsystem, -1)), :]

gdf_ = groupby(_data, :subsystem)
f_ = [SINDy(gdf, [:dV, :dI], [:V, :I], λ = 0.001) for gdf in gdf_]
print(f_[1], ["V", "I"])
print(f_[2], ["V", "I"])

using DecisionTree
Dtree = DecisionTreeClassifier(max_depth=100, n_subfeatures=1)
features = Matrix(__data[:, [:V, :I, :Vr]])
fit!(Dtree, features, __data.subsystem); # print_tree(Dtree, my_depth);
acc = sum(subsystem .== predict(Dtree, Matrix(_data[:, [:V, :I, :Vr]]))) / length(subsystem)

DT = Dtree
v = x
k = 1
predict(DT, [v; anc_[k]])
dt = 1e-7
x = collect(data[1, [:V, :I]])
y = factory_buck(DataFrame, dr.idx, dr.E, x, (0, 0.01))
ŷ = DataFrame(solve(f_, x, dt, y.Vr, Dtree, y.Vr), ["V", "I"]);
plot(xlabel = L"t", ylabel = L"I(t)", size = (800, 300), margin = 5mm)
plot!(y.t[1:100:end], y.I[1:100:end], label = "true", lw = 2)
plot!(y.t[1:100:end], ŷ.I[1:100:end], label = "pred", lw = 2, ls = :dash, color = :red)
png("212")