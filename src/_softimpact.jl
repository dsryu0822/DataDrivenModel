include("../src/DDM.jl")
include("../src/factorio.jl")
include("../src/visual.jl")
using DecisionTree, Random

d_range = 0.1:0.0001:0.3
plan = DataFrame(idx=eachindex(d_range), d=d_range)
data_tspan = [0, 105]

dr = first(eachrow(plan))
@time if "G:/DDM/cached_soft.csv" |> isfile
    data = CSV.read("G:/DDM/cached_soft.csv", DataFrame)
else
    data = factory_soft(DataFrame, dr.idx, dr.d, tspan = data_tspan)
    CSV.write("G:/DDM/cached_soft.csv", data)
end
cutidx = 10_000_000
trng = data[1:cutidx,:]
test = data[cutidx:end,:]

sampled = rand(1:nrow(trng), 15)
f = SINDy(trng[sampled, :], [:dt, :du, :dv], [:t, :u, :v], f_ = [sign, abs, cospi])
print(f, ["t", "u", "v"])

@time error_ = norm.(eachrow(Matrix(trng[:, [:dt, :du, :dv]])) .- f.(eachrow(Matrix(trng[:, [:t, :u, :v]]))))
bit_alien = error_ .> 1e-4
# scatter(log10.(error_)[1:100:100000], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")

subsystem = ones(Int64, nrow(trng));
subsystem[bit_alien] .= 2;
subsystem[1] = 1;
trng[:, :subsystem] = subsystem;
# plot(trng.u[1:100:100000], trng.v[1:100:100000], color=trng.subsystem[1:100:100000], xlabel=L"u", ylabel=L"v", label=:none, ms=1, alpha=0.5, size=(800, 800))

idx_sltd = findall((subsystem .!= circshift(subsystem, 1)) .|| (subsystem .!= circshift(subsystem, -1)))

T_ = trunc.(Int64, exp10.(4:0.1:6))
log_pfmc = DataFrame(zeros(0, 1+length(T_)), :auto)
for sparsity = [1]#, 20, 10, 1]
    println("sparsity = $sparsity")
    len_exact = []
    for T = T_
        print("T = $T: ")
        _trng = trng[(1:sparsity:T) ∪ ((1:T) ∩ idx_sltd),:]
        gdf_ = groupby(_trng, :subsystem)
        f_ = [SINDy(gdf, [:dt, :du, :dv], [:t, :u, :v], f_ = [sign, abs, cospi], λ = 1e-5) for gdf in gdf_]
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
        push!(len_exact, y.t[idx_miss])
        xtk = [100, y.t[idx_miss], 105]
        if T ∈ Int64.(exp10.(1:8))
            plot(xlabel = L"t", ylabel = L"u(t)", size = (800, 200), xlims = [100, 105], xticks = xtk, xformatter = x -> "$(round(x, digits = 5))", margin = 5mm, legend = :none)
            plot!(y.t[1:100:end], y.u[1:100:end], label = "true", lw = 2)
            plot!(y.t[1:100:end], ŷ.u[1:100:end], label = "pred", lw = 2, ls = :dash, color = :red)
            png("T = $T.png")
        end
    end
    push!(log_pfmc, [sparsity; len_exact])
    CSV.write("log_pfmc_soft.csv", log_pfmc)
    scatter(log10.(T_), collect(log_pfmc[end, Not(1)]) .- 1, xlabel = L"\log_{10} T", ylabel = "Perfect predicted", ylims = [0,Inf], legend = :none)
end
cog = range(colorant"orange", colorant"green", length = 6)
plot(xlabel = L"\log_{10} T", ylabel = "break time", ylims = [0,5.1])
plot!(log10.(T_), collect(log_pfmc[1, 2:end]) .- 100, color = cog[1], label = 100/log_pfmc.x1[1], lw = 3, alpha = 0.5)
png("soft.png")