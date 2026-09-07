include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    four wing

''''''''''''''''''''''''''''''''''''''''''''''''''"""

sol = factory_fourwing(DataFrame, 0.135, saveat = 2500:1e-2:3000)
plot(sol.x, sol.y, sol.z, alpha = .5, color = :black, xlabel = L"x", ylabel = L"y", zlabel = L"z", size = [400, 400])

pm, pM = 0.12, 0.16; p0, p1 = 0.135, 0.14;
p_ = range(pm, pM, length = 1001)
bfcn = callbfcn("G:/BF/fourwing/bfcnA.jld2")
@showprogress @threads for k in eachindex(p_)
    sol = factory_fourwing(DataFrame, p_[k], saveat = 2900:1e-2:3000)
    bfcn[p_[k]] = sol.z[arglmax(sol.z)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [pm, p0, p1, pM, .21], ms = .5, ma = .5, msw = 0, color = :black); png("temp")
_bfcn = deepcopy(bfcn)
bfcn = deepcopy(_bfcn)
# JLD2.@save "G:/BF/fourwing/bfcnA.jld2" bfcn

trajA0 = factory_fourwing(DataFrame, p0, saveat = 2900:1e-2:3000)
trajA1 = factory_fourwing(DataFrame, p1, saveat = 2900:1e-2:3000)
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
    sol = solve(ODEProblem(f_, [0.1, 0.1, 0.1], (0, 3000), (β_[k],)), RK4(), dt = 1e-2, adaptive=false, maxiters = Inf)
    z_ = sol[3, sol.t .≥ 2900]
    bfcn[β_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .5, msw = 0, color = :blue); png("temp1")
# JLD2.@save "G:/BF/fourwing/bfcnB.jld2" bfcn
