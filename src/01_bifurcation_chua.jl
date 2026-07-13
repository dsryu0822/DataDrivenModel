include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Chua circuit

''''''''''''''''''''''''''''''''''''''''''''''''''"""
sol = factory_chua(DataFrame, 8.8)
plot(sol.x, sol.y, sol.z, alpha = .5)

pm, pM = 8.430, 8.8; p0, p1 = 8.79, 8.80;
p_ = range(pm, pM, length = 1001)
bfcn = callbfcn("G:/BF/chua/bfcnA.jld2")
@showprogress @threads for k in eachindex(p_)[8.43 .≤ p_ .≤ 8.5]
    sol = factory_chua(DataFrame, p_[k], saveat = 0:1e-3:3000)
    z_ = sol.z[sol.t .≥ 2000]
    bfcn[p_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [pm, p0, p1, pM], ms = .5, ma = .5, msw = 0, color = :black); png("temp")
# JLD2.@save "G:/BF/chua/bfcnA.jld2" bfcn

trajA0 = factory_chua(DataFrame, p0, saveat = 1000:1e-2:1500)
trajA1 = factory_chua(DataFrame, p1, saveat = 1000:1e-2:1500)
vrbl = reverse(half(names(trajA0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:3)
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
CSV.write("G:/BF/chua/trajA0.csv", trajA0)
CSV.write("G:/BF/chua/trajA1.csv", trajA1)
CSV.write("G:/BF/chua/trajB0.csv", trajB0)
CSV.write("G:/BF/chua/trajB1.csv", trajB1)

ic_ = shuffle(trajA0)[1:100, :]
βm = (pm - p0) / (p1 - p0)
βM = (p1 - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 1001)
f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
# bfcn = callbfcn("G:/BF/chua/bfcnB.jld2")
bfcn = callbfcn()
@showprogress @threads for k in eachindex(β_)
    for l in 1:100
        try
            sol = ssolve(f_[k], ic_[[l], f1.rname], 2000:1e-2:3000)
            z_ = sol.z[sol.t .≥ 2000]
            bfcn[β_[k]] = z_[arglmax(z_)]
            break
        catch
            continue
        end
    end
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .5, msw = 0, color = :red); png("temp")
# JLD2.@save "G:/BF/chua/bfcnB.jld2" bfcn

