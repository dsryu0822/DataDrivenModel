include("../src/DDM.jl")
include("../src/factorio.jl")
include("../src/visual.jl")
using DecisionTree, Random

d_range = 0.1:0.0001:0.3
plan = DataFrame(idx=eachindex(d_range), d=d_range)
data_tspan = [0, 105]

dr = first(eachrow(plan))
@time if "C:/DDM/cached_soft.csv" |> isfile
    data = CSV.read("C:/DDM/cached_soft.csv", DataFrame)
else
    data = factory_soft(DataFrame, dr.idx, dr.d, tspan = data_tspan)
    CSV.write("C:/DDM/cached_soft.csv", data)
end

sampled = rand(1:nrow(data), 22)
f = SINDy(data[sampled, :], [:dt, :du, :dv], [:t, :u, :v], M = 3)
print(f, ["t", "u", "v"])

@time error_ = norm.(eachrow(Matrix(data[:, [:dt, :du, :dv]])) .- f.(eachrow(Matrix(data[:, [:t, :u, :v]]))))
bit_alien = error_ .> 1e-4
# scatter(log10.(error_)[1:100:100000], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")

subsystem = ones(Int64, nrow(data));
subsystem[bit_alien] .= 2; subsystem[1] = 1;
data[:, :subsystem] = subsystem;
# plot(data.u[1:100:100000], data.v[1:100:100000], color=data.subsystem[1:100:100000], xlabel=L"u", ylabel=L"v", label=:none, ms=1, alpha=0.5, size=(800, 800))

cutidx = 10_000_000
trng = data[1:cutidx,:]
test = data[cutidx:end,:]

τ_ = trunc.(Int64, exp10.(0.5:0.1:2))
# scatter(findall(trng.subsystem .!= circshift(trng.subsystem, -1))[τ_], τ_, scale = :log10, size = [600,600], xlabel = L"T", ylabel = L"τ", label = :none)
# log_pfmc = DataFrame(zeros(0, length(τ_)), :auto)
len_exact = []
for T = findall(trng.subsystem .!= circshift(trng.subsystem, -1))[τ_]
    print("T = $T: ")

    _trng = trng[1:T,:]
    gdf_ = groupby(_trng, :subsystem)
    f_ = [SINDy(gdf, [:dt, :du, :dv], [:t, :u, :v], M = 5) for gdf in gdf_]
    # print.(f_, Ref(["t", "u", "v"]))

    labels = _trng.subsystem
    features = Matrix(_trng[:, [:t, :u, :v]])
    acc_ = []
    for seed in 1:10
        Dtree = build_tree(_trng.subsystem, features, 2, 2, rng = seed); # print_tree(Dtree, feature_names = ["u", "v"])
        push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break; else print("█") end
    end
    Dtree = build_tree(_trng.subsystem, features, 2, 2, rng = argmax(acc_))
    println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")

    dt = 1e-5
    x = collect(test[1, [:t, :u, :v]])
    y = test
    ŷ = DataFrame(solve(f_, x, dt, y.t, Dtree), ["t", "u", "v"])
    idx_miss = findfirst([abs.(y.u - ŷ.u) .> .1; true]) - 1
    push!(len_exact, count(findall(test.subsystem .!= circshift(test.subsystem, -1)) .< idx_miss))
    xtk = [trunc(y.t[idx_miss], digits = 2)]
    # plot(xlabel = L"t", ylabel = L"u(t)", size = (800, 200), xlims = [100, 105], xticks = xtk, xformatter = x -> "$(round(x, digits = 5))", margin = 5mm, legend = :none)
    plot(legend = :none, size = [600, 200], ticks = false)
    plot!(y.t[1:100:end], y.u[1:100:end], label = "true", lw = 3)
    plot!(y.t[1:100:end], ŷ.u[1:100:end], label = "pred", lw = 3, ls = :dash, color = :red2)
    png("T = $T.png")
end
CSV.write("log_pfmc_soft.csv", DataFrame(input = τ_, output = len_exact))

plot(log10.(τ_), len_exact, legend = :none, lw = 3, alpha = 0.5, ticks = false)
png("soft.png")