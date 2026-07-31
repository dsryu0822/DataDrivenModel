include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    food chain

''''''''''''''''''''''''''''''''''''''''''''''''''"""
sol = factory_foodchain(DataFrame, .99976, ic = [.5, .5, .8], saveat = 0:1e-1:10000)[40000:end,:]
icc = shuffle(sol)

pm, pM = .88, 1.00; p0, p1 = .950, .960;
K_ = range(pm, pM, length = 1001)

trajA0 = factory_foodchain(DataFrame, p0, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
trajA1 = factory_foodchain(DataFrame, p1, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
vrbl = reverse(half(names(trajA0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:4)
f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-8); f0 |> println
f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-8); f1 |> println
trajB0 = ssolve(f0, trajA0[[1], f0.rname], 4000:1e-2:5000)
trajB1 = ssolve(f1, trajA1[[1], f1.rname], 4000:1e-2:5000)
cnfg = cookPI(vrbl, poly = 0:3)
g0 = SINDyPI(trajA0, vrbl, cnfg; λ = 1e-8); g0 |> println
g1 = SINDyPI(trajA1, vrbl, cnfg; λ = 1e-8); g1 |> println
trajC0 = ssolve(g0, trajA0[[1], g0.rname], 4000:1e-2:5000)
trajC1 = ssolve(g1, trajA1[[1], g1.rname], 4000:1e-2:5000)

# JLD2.@load "foodchain_transient.jld2"
ε_ = exp10.(range(-5, -2, 20))
K_c = 0.99976; K_ = K_c .+ ε_
# β_c = 4.612; β_ = β_c .+ 15ε_ # It's possible to recover
β_c = 4.615; β_ = β_c .+ ε_
α_c = 4.777865647; α_ = α_c .+ 60ε_
τA__ = [zeros(100) for K in K_]
τB__ = [zeros(300) for β in β_]
τC__ = [zeros(100) for α in α_]
# JLD2.@save "foodchain_transient_failed.jld2" K_c β_c α_c ε_ τA__ τB__ τC__
# [append!(τA_, zeros(100)) for τA_ in τA__]
# [append!(τB_, zeros(100)) for τB_ in τB__]
# [append!(τC_, zeros(100)) for τC_ in τC__]

for k in eachindex(K_)
    println("K = $(K_[k])")
    @showprogress @threads for l in eachindex(τA__[k])
        if !iszero(τA__[k][l]) continue end
        sol = factory_foodchain(DataFrame, K_[k], ic = [icc[l, [:R, :C, :P]]...], saveat = 0:1e-1:10000; callback = early_stop)
        τ = findlast(sol.P .> 0.55)
        τA__[k][l] = !isnothing(τ) ? sol.t[τ] : last(sol.t)
    end
end
scatter(ε_, mean.(τA__), scale = :log10, ticks = exp10.(-5:6))

cnfg = cook(vrbl, poly = 0:4)
f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
for k in eachindex(β_)
    println("$k: β = $(β_[k])")
    @showprogress @threads for l in eachindex(τB__[k])
        if !iszero(τB__[k][l]) continue end
        sol = ssolve(f_[k], icc[[l], f_[k].rname], 0:1e-1:10000)
        τ = findlast(sol.P .> 0.10)
        τB__[k][l] = !isnothing(τ) ? sol.t[τ] : last(sol.t)
    end
end
scatter(ε_, mean.(τB__), scale = :log10, ticks = exp10.(-5:6))

cnfg = cookPI(vrbl, poly = 0:3)
g_ = [syntheticSINDy((1-α)*g0.matrix + α*g1.matrix, vrbl, cnfg, method = "SINDyPI") for α in α_]
for k in eachindex(α_)
    println("$k: α = $(α_[k])")
    @showprogress @threads for l in eachindex(τC__[k])
        if !iszero(τC__[k][l]) continue end
        sol = ssolve(g_[k], icc[[l], g_[k].rname], 0:1e-1:10000)
        τ = findlast(sol.P .> 0.55)
        τC__[k][l] = !isnothing(τ) ? sol.t[τ] : last(sol.t)
    end
end
scatter(ε_, mean.(τC__), scale = :log10, ticks = exp10.(-5:6))


plot(size = [400, 400], ylims = exp10.([2,4]), legend = true, xlabel = L"\mu - \mu_{c}", ylabel = L"\tau")
scatter!(ε_, mean.(τA__), scale = :log10, ticks = exp10.(-5:6), label = "ground truth", shape = :rect, color = :white)
scatter!(ε_, mean.(τB__), scale = :log10, ticks = exp10.(-5:6), label = "SINDy", shape = :x, color = :red)
scatter!(ε_, mean.(τC__), scale = :log10, ticks = exp10.(-5:6), label = "implicit-SINDy", shape = :+, color = :blue)
png("G:/foodchain_transient.png")

temp = callbfcn("G:/BF/foodchain/bfcnB_9596.jld2")
maximum(keys(temp))
cnfg = cook(vrbl, poly = 0:4)
tempf = syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy")
sol = ssolve(tempf, icc[[7], tempf.rname], 0:1e-1:10000)
plot(sol.R, sol.C, sol.P)
