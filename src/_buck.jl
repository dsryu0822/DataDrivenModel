include("../src/DDM.jl")
include("../src/factorio.jl")
include("../src/visual.jl")
using DecisionTree, Random

E_range = 15:0.001:40
plan = DataFrame(idx=eachindex(E_range), E=E_range)
data_tspan = [0, 1.005]

dr = last(eachrow(plan))
@time if "G:/DDM/cached_buck.csv" |> isfile
    data = CSV.read("G:/DDM/cached_buck.csv", DataFrame)
else
    data = factory_buck(DataFrame, dr.idx, dr.E, tspan = data_tspan)
    CSV.write("G:/DDM/cached_buck.csv", data)
end

sampled = rand(1:nrow(data), 4)
f = SINDy(data[sampled, :], [:dV, :dI], [:V, :I], verbose=true)
print(f, ["V", "I"])

@time error_ = norm.(eachrow(Matrix(data[:, [:dV, :dI]])) .- f.(eachrow(Matrix(data[:, [:V, :I]]))))
bit_alien = error_ .> 1e-4
# scatter(log10.(error_)[1:10000:end], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")

subsystem = ones(Int64, nrow(data));
subsystem[bit_alien] .= 2;
data[:, :subsystem] = subsystem;
# plot(data.V[1:100:100000], data.I[1:100:100000], color=data.subsystem[1:100:100000], xlabel=L"V", ylabel=L"I", label=:none, ms=1, alpha=0.5, size=(800, 800))

cutidx = 10_000_000
trng = data[1:cutidx,:]
test = data[cutidx:end,:]

# idx_sltd = findall((trng.subsystem .!= circshift(trng.subsystem, 1)) .|| (trng.subsystem .!= circshift(trng.subsystem, -1)))

τ_ = unique(trunc.(Int64, exp10.(.5:.1:3)))
# scatter(findall(trng.subsystem .!= circshift(trng.subsystem, -1))[τ_], τ_, scale = :log10, size = [600,600], xlabel = L"T", ylabel = L"τ", label = :none)
# log_pfmc = DataFrame(zeros(0, 1+length(τ_)), :auto)
len_exact = []
for T = findall(trng.subsystem .!= circshift(trng.subsystem, -1))[τ_]
    print("T = $T: ")

    _trng = trng[1:T,:]
    gdf_ = groupby(_trng, :subsystem)
    f_ = [SINDy(gdf, [:dV, :dI], [:V, :I]) for gdf in gdf_]
    # print.(f_, Ref(["V", "I"]))

    labels = _trng.subsystem
    features = Matrix(_trng[:, [:V, :I, :Vr]])
    acc_ = []
    for seed in 1:10
        Dtree = build_tree(_trng.subsystem, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
        push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break; else print("█") end
    end
    Dtree = build_tree(_trng.subsystem, features, rng = argmax(acc_))
    println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")

    dt = 1e-7
    x = collect(test[1, [:V, :I]])
    y = test
    ŷ = DataFrame(solve(f_, x, dt, y.Vr, Dtree, y.Vr), ["V", "I"]);
    idx_miss = findfirst([abs.(y.I - ŷ.I) .> 5e-2; true]) - 1
    push!(len_exact, count(findall(test.subsystem .!= circshift(test.subsystem, -1)) .< idx_miss))
    xtk = [trunc(y.t[idx_miss], digits = 2)]
    plot(legend = :none, size = [600, 200])
    plot!(y.t[1:100:end], y.I[1:100:end], label = "true", lw = 2)
    plot!(y.t[1:100:end], ŷ.I[1:100:end], label = "pred", lw = 2, ls = :dash, color = :red)
    png("T = $T.png")
end
CSV.write("log_pfmc_buck.csv", DataFrame(input = τ_, output = len_exact))

plot(log10.(τ_), len_exact, legend = :none, lw = 3, alpha = 0.5, ticks = false)
png("buck.png")