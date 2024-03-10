# include("../src/factorio.jl")
# include("../src/DDM.jl")
# include("../src/visual.jl")

# data = factory_HR(0, 1, (0,200))
# plot(data.x, data.z, legend = :none, xlabel = "x", ylabel = "z")
# hline!([1, -1], color = :red)
# plot(data.t, data.z, legend = :none, xlabel = "t", ylabel = "dz")
# hline!([1, -1], color = :red)

# foo = vec(sum(abs2, Matrix(data[:, [:x, :y, :z]]), dims = 2))

# plot(abs.((foo ./ circshift(foo, 1)) .- 1)[Not(1)], xlabel = "t", ylabel = "r", title = "algorithm 1")
# bar = abs.((foo ./ circshift(foo, 1)) .- 1)

# plot(data.z, legend = :none, xlabel = "t_k", ylabel = "z", title = "algorithm 1")
# hline!([1, -1], color = :red)
# vline!(findall(bar .> 0.1), alpha = 0.5, color = :blue, style = :dash)

# cd("//155.230.155.221/ty/DS");
# pwd()

# look(data)
# _data = data[1:end, :]; n = nrow(_data)
# # plot(data.V, data.I)

# using ProgressBars

# normeddf = sum.(abs2, eachrow(diff(Matrix(_data[:, [:dx, :dy, :dz]]), dims = 1)))

# # plot(abs.(diff(normeddf)))
# # plot(normeddf)

# # plot(_data.z, xlims = (0, 50000))
# # hline!([1, -1], color = :blue)

# jumpT = [1; findall(normeddf .> 0.02)]
# vline!(jumpT, color = :red)

# T_ = UnitRange.(jumpT, circshift(jumpT, -1)); pop!(T_); T_
# idx_T_ = [argmax(length.(T_))]

# STLSQ_ = []
# for i in axes(T_, 1)
#     sampled = []
#     for j in 1:100
#         sampled = rand(reduce(vcat, T_[[idx_T_; i]]), 200)
#         X = Θ(_data[sampled, [:t, :x, :y, :z]], K = 5)
#         if rank(X) ≥ max(size(X, 2))
#             break
#         end
#     end
#     if rank(X) < max(size(X, 2))
#         println("sampling failed: rank = ", rank(X), " < ", max(size(X, 2)))
#     end
    
#     STLSQed = STLSQ(_data[sampled, :], [:dx, :dy, :dz], [:t, :x, :y, :z], K = 5, λ = 0.001)
#     push!(STLSQ_, STLSQed)
#     println("\ni=$i, MSE=$(STLSQed.MSE) with ", length(STLSQed.matrix.nzval))
# end

# STLSQ_[2]
# argmin(length.(getproperty.(getproperty.(STLSQ_[Not(8)], :matrix), :nzval)))
# Matrix.(getproperty.(STLSQ_[[2, 8]], :matrix))[1]
# Matrix.(getproperty.(STLSQ_[[4, 8]], :matrix))[2]

# Θ(_data[1:10, [:t, :x, :y, :z]], K = 5)

# sampled = rand(reduce(vcat, T_[[idx_T_; 2]]), 200)
# STLSQed = STLSQ(_data[sampled, :], [:dx, :dy, :dz], [:t, :x, :y, :z], K = 5, λ = 0.001)
# α = Matrix(STLSQ_[4].matrix)

# stranger = Int64[]
# error_ = Float64[]
# for k in ProgressBar(1:n)
#     error = sum(abs2, collect(_data[k, [:dx, :dy, :dz]])' .- Θ(_data[k, [:t, :x, :y, :z]], K = 3, f_ = [cos])*α)
#     push!(error_, error)
#     if error > 1e-8
#         push!(stranger, k)
#     end
# end
# scatter(log10.(error_))

# subsystem = ones(Int64, n);
# subsystem[stranger] .= 2;
# _data[:, :subsystem] = subsystem;


# idx_left = findall([count(stranger .∈ Ref(TT)) for TT in T_] .> 100)
# T_ = T_[idx_left]

# idx_T_ = [argmax(length.(T_))]

# STLSQ_ = []
# for i in axes(T_, 1)
#     sampled = []
#     for j in 1:100
#         sampled = rand(reduce(vcat, T_[[idx_T_; i]]), 50)
#         if rank(Θ(_data[sampled, [:t, :x, :y, :z]], K = 3, f_ = [cos])) ≥ 39
#             break
#         end
#     end
    
#     STLSQed = STLSQ(_data[sampled, :], [:dx, :dy, :dz], [:t, :x, :y, :z], K = 3, f_ = [cos], λ = 0.001)
#     push!(STLSQ_, STLSQed)
#     println("\ni=$i, MSE=$(STLSQed.MSE) with ", length(STLSQed.matrix.nzval))
# end

# α = Matrix(STLSQ_[3].matrix)

# stranger = Int64[]
# error_ = Float64[]
# for k in ProgressBar(1:n)
#     error = sum(abs2, collect(_data[k, [:dx, :dy, :dz]])' .- Θ(_data[k, [:t, :x, :y, :z]], K = 3, f_ = [cos])*α)
#     push!(error_, error)
#     if error > 1e-8
#         push!(stranger, k)
#     end
# end
# scatter(log10.(error_))
# _data.subsystem[findall(subsystem .> 1) ∩ stranger] .= 3;


# plot(_data.x[1:10:end], _data.z[1:10:end],
#     xlabel=L"x", ylabel=L"z",
#     color=_data.subsystem[1:10:end], label=:none,
#     size=(800, 800))
# hline!([1, -1], color = :black, style = :dash)
# png("temp 1")

# gdf_ = groupby(_data, :subsystem)
# STLSQ_ = [STLSQ(gdf, [:dx, :dy, :dz], [:t, :x, :y, :z], K = 3, f_ = [cos]) for gdf in gdf_]
# Matrix.(getproperty.(STLSQ_, :matrix))[1]
# Matrix.(getproperty.(STLSQ_, :matrix))[2]
# Matrix.(getproperty.(STLSQ_, :matrix))[3]

# function factory_STLSQ(STLSQed)
#     function f(s, x)
#         return [1; vec(Θ(x, K = 3, f_ = [cos]) * STLSQed[s].matrix)]
#     end
#     return f
# end
# g = factory_STLSQ(STLSQ_)


# using DecisionTree
# my_depth = 3
# Dtree = DecisionTreeClassifier(max_depth=my_depth)
# features = Matrix(_data[:, [:t, :x, :y, :z]])
# fit!(Dtree, features, _data.subsystem)
# print_tree(Dtree, my_depth)
# acc = sum(_data.subsystem .== predict(Dtree, features)) / length(subsystem)
# println(acc)

# x_ = [collect(data[1, [:t, :x, :y, :z]])]
# x = x_[end]

# dt = 10^(-3)
# for t in ProgressBar(0:dt:200)
#     s = predict(Dtree, x_[end])
#     x, dx = RK4(g, s, x_[end], dt)
#     push!(x_, x)
# end
# _x_ = stack(x_)[:, 1:(end-2)]

# uv1 = plot(data.t[1:10:end], data.x[1:10:end], ylabel = L"x", label="data")
# plot!(uv1, data.t[1:10:end], _x_[2, 1:10:end], color=:red, alpha = 0.5, label="predicted")
# title!(uv1, "Hindmarsh-Rose neuronal model")

# uv2 = plot(data.t[1:10:end], data.z[1:10:end], ylabel = L"z", label="data")
# plot!(uv2, data.t[1:10:end], _x_[4, 1:10:end], color=:red, alpha= 0.5, label="predicted")
# xlabel!(uv2, L"t")
# plot(uv1, uv2, layout=(2, 1), size=(800, 800));
# png("temp 2");

# # ---

# @time data = factory_buck(dr.idx, dr.E, (0.495, 0.55))
# x_ = [collect(data[1, [:V, :I]])]
# x = x_[end]

# for ramp in ProgressBar(data.Vr)
#     s = predict(Dtree, [x_[end]; ramp])
#     x, dx = RK4(g, s, x_[end], dt)
#     push!(x_, x)
# end
# _x_ = stack(x_)

# uv1 = plot(data.t[1:10:end], data.V[1:10:end], ylabel = L"V", label="data")
# plot!(uv1, data.t[1:10:end], _x_[1, 1:10:end], color=:red, style=:dash, label="predicted")
# title!(uv1, "Buck converter")

# uv2 = plot(data.t[1:10:end], data.I[1:10:end], ylabel = L"I", label="data")
# plot!(uv2, data.t[1:10:end], _x_[2, 1:10:end], color=:red, style=:dash, label="predicted")
# xlabel!(uv2, L"t")
# plot(uv1, uv2, layout=(2, 1), size=(1600, 900))


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
trng = data[1:100000,:]
trng = data[1:cutidx,:]
test = data[cutidx:end,:]

normeddf = sum.(abs2, eachrow(diff(Matrix(trng[:, [:dx, :dy, :dz]]), dims = 1)))

jumpt = [1; findall(normeddf .> 0.02)]
sets = collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1))); pop!(sets)
# plot(trng.t, trng.z)
# hline!([1, -1], color = :blue)
# vline!(trng.t[jumpt], color = :red, ls = :dash)

subsystem = zeros(Int64, nrow(trng));

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
bit_alien = error_ .> 1e-4
# scatter(log10.(error_)[1:100:100000], ylabel = L"\log_{10} | r |", title = "Residuals", legend = :none, xlabel = "Index")
subsystem[.!bit_alien] .= 1;
# plot(trng.u[1:100:100000], trng.v[1:100:100000], color=trng.subsystem[1:100:100000], xlabel=L"u", ylabel=L"v", label=:none, ms=1, alpha=0.5, size=(800, 800))

sets = filter(!isempty, [setdiff(set, findall(.!bit_alien)) for set in sets])
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
bit_alien = error_ .> 1e-4
subsystem[bit_alien] .= 2;

trng[:, :subsystem] = subsystem;


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