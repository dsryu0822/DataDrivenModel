include.("../core/" .* readdir("core")[[1,2,3,4,6]])

trajA0 = factory_rossler(DataFrame, 14, saveat = 900:1e-3:1000)
trajA1 = factory_rossler(DataFrame, 15, saveat = 900:1e-3:1000)
plot(plot(trajA0.x, trajA0.y), plot(trajA1.x, trajA1.y))

vrbl = (["dx", "dy"], ["x", "y"])
cnfg = cook(vrbl, poly = 0:4)
f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-8); f0 |> println
f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-8); f1 |> println
trajB0 = ssolve(f0, trajA0[[1], f0.rname], 900:1e-3:1000)
trajB1 = ssolve(f1, trajA1[[1], f1.rname], 900:1e-3:1000)

plot(plot(trajB0.x, trajB0.y), plot(trajB1.x, trajB1.y))

