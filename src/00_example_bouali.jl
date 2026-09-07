include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Bouali

''''''''''''''''''''''''''''''''''''''''''''''''''"""

sol = factory_bouali(DataFrame, 0.385, ic = [2,1,0], saveat = 4900:1e-2:5000)
plot(sol.x, sol.y, sol.z, alpha = .5, color = :black, xlabel = L"x", ylabel = L"y", zlabel = L"z", size = [400, 400])

pm, pM = 0.5, 1.2; p0, p1 = 0.7, 0.8;
p_ = range(pm, pM, length = 1001)
bfcn = callbfcn("G:/BF/bouali/bfcnA.jld2")
@showprogress @threads for k in eachindex(p_)
    sol = factory_bouali(DataFrame, p_[k], saveat = 2500:1e-2:3000)
    bfcn[p_[k]] = sol.x[arglmax(sol.x)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [pm, p0, p1, pM], ms = .5, ma = .5, msw = 0, color = :red); png("temp1")
scatter(dict2bifurcation(bfcn)..., xlims = [0.378, 0.4], ms = .5, ma = .5, msw = 0, color = :red); png("temp1")
# JLD2.@save "G:/BF/bouali/bfcnA.jld2" bfcn

trajA0 = factory_bouali(DataFrame, p0, saveat = 2500:1e-2:3000)
trajA1 = factory_bouali(DataFrame, p1, saveat = 2500:1e-2:3000)
vrbl = reverse(half(names(trajA0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:4)
f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f0 |> println
f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f1 |> println
trajB0 = ssolve(f0, trajA0[[1], f0.rname], 1000:1e-2:1500)
trajB1 = ssolve(f1, trajA1[1, f1.rname], 1000:1e-2:1500)
plot(
    plot(trajA0.x, trajA0.y, trajA0.z, alpha = .5, color = :black),
    plot(trajA1.x, trajA1.y, trajA1.z, alpha = .5, color = :black),
    plot(trajB0.x, trajB0.y, trajB0.z, alpha = .5, color = :red),
    plot(trajB1.x, trajB1.y, trajB1.z, alpha = .5, color = :red),
)

βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 1001)
f_ = affine(Function, f0, f1) # affine(String, f0, f1) |> println
bfcn = callbfcn()
@showprogress @threads for k in eachindex(β_)
    sol = solve(ODEProblem(f_, [1.0, 1.0, 0.1], (0, 3000), (β_[k],)), RK4(), dt = 1e-2, adaptive=false, maxiters = Inf)
    x_ = sol[1, sol.t .≥ 2500]
    bfcn[β_[k]] = x_[arglmax(x_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .5, msw = 0, color = :blue); png("temp1")
# JLD2.@save "G:/BF/bouali/bfcnB.jld2" bfcn
