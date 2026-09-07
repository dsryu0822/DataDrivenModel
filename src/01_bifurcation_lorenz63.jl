include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Lorenz-63

''''''''''''''''''''''''''''''''''''''''''''''''''"""
sol = factory_lorenz63(DataFrame, [10, 28, 8/3])
plot(sol.x, sol.y, sol.z, alpha = .5)

pm, pM = -2, 3; p0, p1 = 0, 1;
p_ = range(pm, pM, length = 2001)
σ_ = range(6, 15, length = 2001)
ρ_ = range(120, 150, length = 2001)
b_ = range(3, 5, length = 2001)
bfcn = callbfcn()
@showprogress @threads for k in eachindex(p_)
    sol = factory_lorenz63(DataFrame, [σ_[k], ρ_[k], b_[k]], ic = [100, 100, 100], saveat = 900:1e-3:1000)
    z_ = sol.z[sol.t .≥ 900]
    bfcn[p_[k]] = z_[arglmax(z_)]
end
# bfcn = callbfcn("G:/BF/lorenz63/bfcnA.jld2")
scatter(dict2bifurcation(bfcn)..., xlims = [pm, pM], xticks = [pm, p0, p1, pM], ms = .1, ma = 1.0, msw = 0, color = :black, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> ""); png("temp")
scatter([zeros(344); ones(381)], [bfcn[0]; bfcn[1]], xlims = [pm, pM], xticks = [pm, p0, p1, pM], ms = 3, ma = .5, msw = 0, color = :black, ylims = [160, 260], yticks = [160, 260], size = [600, 100], shape = :x, xformatter = _ -> ""); png("temp")
# JLD2.@save "G:/BF/lorenz63/bfcnA.jld2" bfcn

trajA0 = factory_lorenz63(DataFrame, [σ_[ 801], ρ_[ 801], b_[ 801]], saveat = 900:1e-3:1000)
trajA1 = factory_lorenz63(DataFrame, [σ_[1201], ρ_[1201], b_[1201]], saveat = 900:1e-3:1000)
vrbl = reverse(half(names(trajA0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:2)
f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f0 |> println
f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f1 |> println
trajB0 = ssolve(f0, trajA0[[1], f0.rname], 900:1e-3:1000)
trajB1 = ssolve(f1, trajA1[[1], f1.rname], 900:1e-3:1000)
cnfg = cookPI(vrbl, poly = 0:2)
g0 = SINDyPI(trajA0, vrbl, cnfg; λ = 1e-3); g0 |> println
g1 = SINDyPI(trajA1, vrbl, cnfg; λ = 1e-3); g1 |> println
trajC0 = ssolve(g0, trajA0[[1], g0.rname], 900:1e-3:1000)
trajC1 = ssolve(g1, trajA1[[1], g1.rname], 900:1e-3:1000)

plot(
    plot(trajA0.x, trajA0.y, trajA0.z, alpha = .5, color = :black),
    plot(trajA1.x, trajA1.y, trajA1.z, alpha = .5, color = :black),
    plot(trajB0.x, trajB0.y, trajB0.z, alpha = .5, color = :red),
    plot(trajB1.x, trajB1.y, trajB1.z, alpha = .5, color = :red),
    plot(trajC0.x, trajC0.y, trajC0.z, alpha = .5, color = :blue),
    plot(trajC1.x, trajC1.y, trajC1.z, alpha = .5, color = :blue),
    ticks = [], layout = (:, 2), size = (400, 600)
)

plot(
    plot(trajA0.x, trajA0.y, trajA0.z, alpha = .5, color = :black, xticks = [minimum(trajA0.x)-2.5], yticks = [maximum(trajA0.y)+4.5], zticks = [minimum(trajA0.z)-3],),
    plot(trajA1.x, trajA1.y, trajA1.z, alpha = .5, color = :black, xticks = [minimum(trajA1.x)-2.5], yticks = [maximum(trajA1.y)+4.5], zticks = [minimum(trajA1.z)-3],),
    formatter = _ -> "", layout = (:, 2), size = (400, 200)
); png("temp")
plot(
    plot(trajB0.x, trajB0.y, trajB0.z, alpha = .5, color = :blue, xticks = [minimum(trajB0.x)-2.5], yticks = [maximum(trajB0.y)+4.5], zticks = [minimum(trajB0.z)-3],),
    plot(trajB1.x, trajB1.y, trajB1.z, alpha = .5, color = :blue, xticks = [minimum(trajB1.x)-2.5], yticks = [maximum(trajB1.y)+4.5], zticks = [minimum(trajB1.z)-3],),
    formatter = _ -> "", layout = (:, 2), size = (400, 200)
); png("temp")

CSV.write("G:/BF/lorenz63/trajA0.csv", trajA0)
CSV.write("G:/BF/lorenz63/trajA1.csv", trajA1)
CSV.write("G:/BF/lorenz63/trajB0.csv", trajB0)
CSV.write("G:/BF/lorenz63/trajB1.csv", trajB1)
CSV.write("G:/BF/lorenz63/trajC0.csv", trajC0)
CSV.write("G:/BF/lorenz63/trajC1.csv", trajC1)

βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 2001)
f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
bfcn = callbfcn()
@showprogress for k in eachindex(β_)
    sol = ssolve(f_[k], trajA0[[1], f_[k].rname], 900:1e-3:1000)
    z_ = sol.z[sol.t .≥ 900]
    bfcn[β_[k]] = z_[arglmax(z_)]
end
# bfcn = callbfcn("G:/BF/lorenz63/bfcnB.jld2")

scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .1, ma = 1.0, msw = 0, color = :blue, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> ""); png("temp")
# JLD2.@save "G:/BF/lorenz63/bfcnB.jld2" bfcn

g_ = [syntheticSINDy((1-β)*g0.matrix + β*g1.matrix, vrbl, cnfg, method = "SINDyPI") for β in β_]
bfcn = callbfcn("G:/BF/lorenz63/bfcnC.jld2")
@showprogress @threads for k in eachindex(β_)
    sol = ssolve(g_[k], trajA0[[1], g_[k].rname], 900:1e-3:1000)
    z_ = sol.z[sol.t .≥ 900]
    bfcn[β_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .5, msw = 0, color = :red, ylims = [160, 260], yticks = [160, 260], size = [400, 100], xformatter = _ -> ""); png("temp")
# JLD2.@save "G:/BF/lorenz63/bfcnC.jld2" bfcn


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    time-dependent

''''''''''''''''''''''''''''''''''''''''''''''''''"""

affine(f0, f1) |> println

f_t = string2function(replace(affine(f0, f1), "β = param[1]" => "β = 1e-3tau"))
sol = solve(ODEProblem(f_t, [1, 1, 1], (-2000, 3000)), RK4(), dt = 1e-3, adaptive=false, maxiters = Inf)

scatter(sol.t[arglmax(sol[3, :])], sol[3, arglmax(sol[3, :])], msw = 0, ms = 1)
png("G:/lorenz63 time dependent.png")