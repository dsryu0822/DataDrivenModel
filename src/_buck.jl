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
cutidx = 10_000_000
trng = data[1:cutidx,:]
test = data[cutidx:end,:]

sampled = rand(1:nrow(trng), 4)
f = SINDy(trng[sampled, :], [:dV, :dI], [:V, :I], verbose=true)
print(f, ["V", "I"])

@time error_ = norm.(eachrow(Matrix(trng[:, [:dV, :dI]])) .- f.(eachrow(Matrix(trng[:, [:V, :I]]))))
bit_alien = error_ .> 1e-4
# scatter(log10.(error_)[1:10000:end], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")

subsystem = ones(Int64, nrow(trng));
subsystem[bit_alien] .= 2;
trng[:, :subsystem] = subsystem;
# plot(trng.V[1:100:100000], trng.I[1:100:100000], color=trng.subsystem[1:100:100000], xlabel=L"V", ylabel=L"I", label=:none, ms=1, alpha=0.5, size=(800, 800))

idx_chgd = findall((subsystem .!= circshift(subsystem, 1)) .|| (subsystem .!= circshift(subsystem, -1)))

T_ = trunc.(Int64, exp10.(4:.5:7))
log_pfmc = DataFrame(zeros(0, 1+length(T_)), :auto)
for sparsity = [1000, 200, 100, 20, 10, 1]
    println("sparsity = $sparsity")
    len_exact = []
    for T = T_           
        print("T = $T: ")

        _trng = trng[(1:sparsity:T) ∪ ((1:T) ∩ idx_sltd),:]
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
        idx_miss = findfirst([(abs.(y.I - ŷ.I) ./ abs.(y.I)) .> .1; true]) - 1
        push!(len_exact, y.t[idx_miss])
        xtk = [1.0, y.t[idx_miss], 1.005]
        if T ∈ Int64.(exp10.(1:8))
            plot(xlabel = L"t", ylabel = L"I(t)", size = (800, 200), xlims = [1.0, 1.005], xticks = xtk, xformatter = x -> "$(round(x, digits = 5))", margin = 5mm, legend = :none)
            plot!(y.t[1:100:end], y.I[1:100:end], label = "true", lw = 2)
            plot!(y.t[1:100:end], ŷ.I[1:100:end], label = "pred", lw = 2, ls = :dash, color = :red)
            png("T = $T.png")
        end
    end
    push!(log_pfmc, [sparsity; len_exact])
    CSV.write("log_pfmc_buck.csv", log_pfmc)
    scatter(log10.(T_), collect(log_pfmc[end, Not(1)]) .- 1, xlabel = L"\log_{10} T", ylabel = "Perfect predicted", ylims = [0,Inf], legend = :none)
end

cog = range(colorant"orange", colorant"green", length = 6)
plot(xlabel = L"\log_{10} T", ylabel = "break time", ylims = [0,0.0051])
plot!(log10.(T_), collect(log_pfmc[6, 2:end]) .- 1, color = cog[6], label = 100/log_pfmc.x1[6], lw = 3, alpha = 0.5)
plot!(log10.(T_), collect(log_pfmc[5, 2:end]) .- 1, color = cog[5], label = 100/log_pfmc.x1[5], lw = 3, alpha = 0.5)
plot!(log10.(T_), collect(log_pfmc[4, 2:end]) .- 1, color = cog[4], label = 100/log_pfmc.x1[4], lw = 3, alpha = 0.5)
plot!(log10.(T_), collect(log_pfmc[3, 2:end]) .- 1, color = cog[3], label = 100/log_pfmc.x1[3], lw = 3, alpha = 0.5)
plot!(log10.(T_), collect(log_pfmc[2, 2:end]) .- 1, color = cog[2], label = 100/log_pfmc.x1[2], lw = 3, alpha = 0.5)
plot!(log10.(T_), collect(log_pfmc[1, 2:end]) .- 1, color = cog[1], label = 100/log_pfmc.x1[1], lw = 3, alpha = 0.5)
png("buck.png")
# scatter!(fill(7, 6), log_pfmc[:, end] .- 1, text = ["6", "5", "    4", "3", "2", "    1"], color = cog, alpha = 0.5, ms = 10, label = :none)
