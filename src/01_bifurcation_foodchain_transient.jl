include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    food chain

''''''''''''''''''''''''''''''''''''''''''''''''''"""
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

ε_ = exp10.(range(-5, -2, 20))

K_c = 0.99976
K_ = K_c .+ ε_

β_c = 0
β_ = β_c .+ ε_
f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]

α_c = 4.777865647
α_ = α_c .+ ε_
g_ = [syntheticSINDyPI((1-α)*g0.matrix + α*g1.matrix, vrbl, cnfg, method = "SINDyPI") for α in α_]
