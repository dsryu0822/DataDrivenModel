include("../src/DDM.jl")
include("../src/factorio.jl")
include("../src/visual.jl")
using DecisionTree, Random

dt = 1e-3
data_tspan = [0, 10050]

@time if "C:/DDM/cached_hrnm.csv" |> isfile
    data = CSV.read("C:/DDM/cached_hrnm.csv", DataFrame)
else
    data = factory_HR(DataFrame, 0, 1, tspan = data_tspan)
    CSV.write("C:/DDM/cached_hrnm.csv", data)
end
cutidx = 10_000_000
trng = data[1:1_000_000,:]
trng = data[1:cutidx,:]
test = data[cutidx:end,:]

normeddf = sum.(abs2, eachrow(diff(Matrix(trng[:, [:dx, :dy, :dz]]), dims = 1)))
jumpt = [1; findall(normeddf .> 0.02)]
# plot(trng.t, trng.z)
# hline!([1, -1], color = :blue)
# vline!(trng.t[jumpt], color = :red, ls = :dash)

subsystem = zeros(Int64, nrow(trng));
sets = collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1))); pop!(sets);
for id_subsys in 1:4
    idx_long = argmax(length.(sets))

    candy_ = []
    head = trng[sets[idx_long], [:t, :x, :y, :z]]
    for i in eachindex(sets)
        body = trng[sets[i], [:t, :x, :y, :z]]
        X = Θ([head; body], N = 3, f_ = [cos])
        if rank(X) ≥ max(size(X, 2))
            candy = SINDy(trng[reduce(vcat, sets[[idx_long, i]]), :],
                    [:dt, :dx, :dy, :dz], [:t, :x, :y, :z],
                    N = 3, f_ = [cos])
            push!(candy_, i => candy.MSE)
            if candy.MSE < 1e-28 break end
        end
    end
    picked = first(candy_[argmin(last.(candy_))])
    f = SINDy(trng[reduce(vcat, sets[[idx_long; picked]]), :],
        [:dt, :dx, :dy, :dz], [:t, :x, :y, :z],
        N = 3, f_ = [cos])
    print(f, ["t", "x", "y", "z"])

    @time error_ = norm.(eachrow(Matrix(trng[:, [:dt, :dx, :dy, :dz]])) .- f.(eachrow(Matrix(trng[:, [:t, :x, :y, :z]]))))
    bit_alien = error_ .> 1e-8
    # scatter(log10.(error_)[1:100:100000], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")
    idx_alien = findall(.!bit_alien)
    subsystem[idx_alien] .= id_subsys
    # plot(trng.u[1:100:100000], trng.v[1:100:100000], color=trng.subsystem[1:100:100000], xlabel=L"u", ylabel=L"v", label=:none, ms=1, alpha=0.5, size=(800, 800))

    # sets = filter(!isempty, [setdiff(set, idx_alien) for set in sets])
    @time sets = filter(x -> isdisjoint(x, idx_alien), sets)
    if isempty(sets)
        @info "all data points are exhausted!"
        break
    end
end
trng[:, :subsystem] = subsystem;

idx_sltd = findall((trng.subsystem .!= circshift(trng.subsystem, 1)) .|| (trng.subsystem .!= circshift(trng.subsystem, -1)))


T_ = trunc.(Int64, exp10.(4:1:6))
log_pfmc = DataFrame(zeros(0, 1+length(T_)), :auto)
for sparsity = [500]#, 20, 10, 1]
    println("sparsity = $sparsity")
    len_exact = []
    for T = T_
# sparsity = 1
# T = 100000
        print("T = $T: ")
        _trng = trng[(1:sparsity:T) ∪ ((1:T) ∩ idx_sltd),:]
        gdf_ = groupby(_trng, :subsystem)
        f_ = [SINDy(gdf, [:dt, :dx, :dy, :dz], [:t, :x, :y, :z], N = 3, f_ = [cos]) for gdf in gdf_]
        print.(f_, Ref(["t", "x", "y", "z"]))

        labels = _trng.subsystem
        features = Matrix(_trng[:, [:t, :x, :y, :z]])
        acc_ = []
        for seed in 1:10
            Dtree = build_tree(_trng.subsystem, features, 2, 2, rng = seed); # print_tree(Dtree, feature_names = ["u", "v"])
            push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
            if maximum(acc_) ≈ 1 break; else print("█") end
        end
        Dtree = build_tree(_trng.subsystem, features, 2, 2, rng = argmax(acc_))
        println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")

        dt = 1e-3
        x = collect(test[1, [:t, :x, :y, :z]])
        y = test
        ŷ = DataFrame(solve(f_, x, dt, y.t, Dtree), ["t", "x", "y", "z"])
        idx_miss = findfirst([abs.(y.y - ŷ.y) .> .1; true]) - 1
        push!(len_exact, y.t[idx_miss])
        xtk = [10000, y.t[idx_miss], 10050]
        if T ∈ Int64.(exp10.(1:8))
            plot(xlabel = L"t", ylabel = L"u(t)", size = (800, 200), xlims = [10000, 10050], xticks = xtk, xformatter = x -> "$(round(x, digits = 5))", margin = 5mm, legend = :none)
            plot!(y.t[1:100:end], y.y[1:100:end], label = "true", lw = 2)
            plot!(y.t[1:100:end], ŷ.y[1:100:end], label = "pred", lw = 2, ls = :dash, color = :red)
            png("T = $T.png")
        end
    end
    push!(log_pfmc, [sparsity; len_exact])
    CSV.write("log_pfmc_hrnm.csv", log_pfmc)
    scatter(log10.(T_), collect(log_pfmc[end, Not(1)]) .- 1, xlabel = L"\log_{10} T", ylabel = "Perfect predicted", ylims = [0,Inf], legend = :none)
end
cog = range(colorant"orange", colorant"green", length = 2)
plot(xlabel = L"\log_{10} T", ylabel = "break time", ylims = [0,50.1])
plot!(log10.(T_), collect(log_pfmc[1, 2:end]) .- 10000, color = cog[1], label = 100/log_pfmc.x1[1], lw = 3, alpha = 0.5)
plot!(log10.(T_), collect(log_pfmc[2, 2:end]) .- 10000, color = cog[2], label = 100/log_pfmc.x1[2], lw = 3, alpha = 0.5)
plot!(log10.(T_), collect(log_pfmc[3, 2:end]) .- 10000, color = cog[3], label = 100/log_pfmc.x1[3], lw = 3, alpha = 0.5)
png("hrnm.png")