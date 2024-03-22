include("../src/DDM.jl")
include("../src/factorio.jl")
include("../src/visual.jl")
using DecisionTree, Random

data_tspan = [0, 10050]

@time if "C:/DDM/cached_hrnm.csv" |> isfile
    data = CSV.read("C:/DDM/cached_hrnm.csv", DataFrame)
else
    data = factory_HR(DataFrame, 0, 1, tspan = data_tspan)
    CSV.write("C:/DDM/cached_hrnm.csv", data)
end

normeddf = sum.(abs2, eachrow(diff(Matrix(data[:, [:dx, :dy, :dz]]), dims = 1)))
jumpt = [1; findall(normeddf .> 0.02)]
# plot(trng.t, trng.z)
# hline!([1, -1], color = :blue)
# vline!(trng.t[jumpt], color = :red, ls = :dash)

subsystem = zeros(Int64, nrow(data));
sets = collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1))); pop!(sets); sets = UnitRange.(first.(sets), last.(sets))
for id_subsys in 1:4
    idx_long = argmax(length.(sets))

    candy_ = []
    head = data[sets[idx_long], [:t, :x, :y, :z]]
    @time for i in eachindex(sets)
        body = data[sets[i], [:t, :x, :y, :z]]
        candy = SINDy(data[reduce(vcat, sets[[idx_long, i]]), :],
                [:dt, :dx, :dy, :dz], [:t, :x, :y, :z],
                N = 3, f_ = [cos])
        X = Θ([head; body], N = 3, f_ = [cos])
        if rank(X) ≥ max(size(X, 2))
            push!(candy_, i => candy.MSE)
            if candy.MSE < 1e-28 break end
        end
    end
    picked = first(candy_[argmin(last.(candy_))])
    f = SINDy(data[reduce(vcat, sets[[idx_long; picked]]), :],
        [:dt, :dx, :dy, :dz], [:t, :x, :y, :z],
        N = 3, f_ = [cos])
    print(f, ["t", "x", "y", "z"])

    @time error_ = norm.(eachrow(Matrix(data[:, [:dt, :dx, :dy, :dz]])) .- f.(eachrow(Matrix(data[:, [:t, :x, :y, :z]]))))
    bit_alien = error_ .> 1e-8
    # scatter(log10.(error_)[1:100:100000], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")
    idx_alien = findall(.!bit_alien)
    subsystem[idx_alien] .= id_subsys
    # plot(data.u[1:100:100000], data.v[1:100:100000], color=data.subsystem[1:100:100000], xlabel=L"u", ylabel=L"v", label=:none, ms=1, alpha=0.5, size=(800, 800))

    # sets = filter(!isempty, [setdiff(set, idx_alien) for set in sets])
    # sets = filter(x -> isdisjoint(x, idx_alien), sets)
    idx_unlabled = []
    for i in eachindex(sets)
        if first(sets[i]) ∉ idx_alien
            push!(idx_unlabled, i)
        end
    end
    sets = sets[idx_unlabled]
    if isempty(sets)
        @info "all data points are exhausted!"
        break
    end
end
data[:, :subsystem] = subsystem;

cutidx = 1_000_000
trng = data[1:cutidx,:]
test = data[cutidx:end,:]

τ_ = trunc.(Int64, exp10.(0.5:0.5:3))
len_exact = []
for T = τ_
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
    xtk = [y.t[idx_miss]]
        plot(xlabel = L"t", ylabel = L"u(t)", size = (800, 200), xlims = [10000, 10050], xticks = xtk, xformatter = x -> "$(round(x, digits = 5))", margin = 5mm, legend = :none)
        plot!(y.t[1:100:end], y.y[1:100:end], label = "true", lw = 2)
        plot!(y.t[1:100:end], ŷ.y[1:100:end], label = "pred", lw = 2, ls = :dash, color = :red)
        png("T = $T.png")
end
CSV.write("log_pfmc_hrnm.csv", DataFrame(input = τ_, output = len_exact))

plot(log10.(τ_), len_exact, legend = :none, lw = 3, alpha = 0.5, ticks = false)
png("hrnm.png")