include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Lorenz-63

''''''''''''''''''''''''''''''''''''''''''''''''''"""
sol = factory_lorenz63(DataFrame, [10, 28, 8/3])
plot(sol.x, sol.y, sol.z, alpha = .5)

pm, pM = 0, 10; p0, p1 = 0, 1;
p_ = range(pm, pM, length = 1001)
σ_ = range(6, 15, length = 1001)
ρ_ = range(120, 150, length = 1001)
b_ = range(3, 5, length = 1001)
bfcn = callbfcn("G:/BF/lorenz63/bfcnA.jld2")
@showprogress @threads for k in eachindex(p_)
    sol = factory_lorenz63(DataFrame, [σ_[k], ρ_[k], b_[k]], ic = [100, 100, 100], saveat = 900:1e-3:1000)
    z_ = sol.z[sol.t .≥ 900]
    bfcn[p_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [pm, p0, p1, pM], ms = .5, ma = .5, msw = 0, color = :black); png("temp")
# JLD2.@save "G:/BF/lorenz63/bfcnA.jld2" bfcn

trajA0 = factory_lorenz63(DataFrame, [σ_[1], ρ_[1], b_[1]], saveat = 900:1e-3:1000)
trajA1 = factory_lorenz63(DataFrame, [σ_[101], ρ_[101], b_[101]], saveat = 900:1e-3:1000)
vrbl = reverse(half(names(trajA0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:2)
f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f0 |> println
f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f1 |> println
trajB0 = ssolve(f0, trajA0[[1], f0.rname], 1000:1e-2:1500)
trajB1 = ssolve(f1, trajA1[[1], f1.rname], 1000:1e-2:1500)
cnfg = cookPI(vrbl, poly = 0:2)
g0 = SINDyPI(trajA0, vrbl, cnfg; λ = 1e-3); g0 |> println
g1 = SINDyPI(trajA1, vrbl, cnfg; λ = 1e-3); g1 |> println
trajC0 = ssolve(g0, trajA0[[1], g0.rname], 1000:1e-2:1500)
trajC1 = ssolve(g1, trajA1[[1], g1.rname], 1000:1e-2:1500)

plot(
    plot(trajA0.x, trajA0.y, trajA0.z, alpha = .5, color = :black),
    plot(trajA1.x, trajA1.y, trajA1.z, alpha = .5, color = :black),
    plot(trajB0.x, trajB0.y, trajB0.z, alpha = .5, color = :red),
    plot(trajB1.x, trajB1.y, trajB1.z, alpha = .5, color = :red),
    plot(trajC0.x, trajC0.y, trajC0.z, alpha = .5, color = :blue),
    plot(trajC1.x, trajC1.y, trajC1.z, alpha = .5, color = :blue)
)
CSV.write("G:/BF/lorenz63/trajA0.csv", trajA0)
CSV.write("G:/BF/lorenz63/trajA1.csv", trajA1)
CSV.write("G:/BF/lorenz63/trajB0.csv", trajB0)
CSV.write("G:/BF/lorenz63/trajB1.csv", trajB1)
CSV.write("G:/BF/lorenz63/trajC0.csv", trajC0)
CSV.write("G:/BF/lorenz63/trajC1.csv", trajC1)

βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 1001)

f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
bfcn = callbfcn("G:/BF/lorenz63/bfcnB.jld2")
@showprogress for k in eachindex(β_)
    sol = ssolve(f_[k], trajA0[[1], f_[k].rname], 2000:1e-4:3000)
    z_ = sol.z[sol.t .≥ 2000]
    bfcn[β_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], ms = .5, ma = .5, msw = 0, color = :red); png("temp")
# JLD2.@save "G:/BF/lorenz63/bfcnB.jld2" bfcn

g_ = [syntheticSINDy((1-β)*g0.matrix + β*g1.matrix, vrbl, cnfg, method = "SINDyPI") for β in β_]
bfcn = callbfcn("G:/BF/lorenz63/bfcnC.jld2")
@showprogress @threads for k in eachindex(β_)
    sol = ssolve(g_[k], trajA0[[1], g_[k].rname], 900:1e-3:1000)
    z_ = sol.z[sol.t .≥ 900]
    bfcn[β_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .5, msw = 0, color = :red); png("temp")
# JLD2.@save "G:/BF/lorenz63/bfcnC.jld2" bfcn

