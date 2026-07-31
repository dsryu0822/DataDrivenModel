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
scatter([zeros(344); ones(381)], [bfcn[0]; bfcn[1]], xlims = [pm, pM], xticks = [pm, p0, p1, pM], ms = 3, ma = .5, msw = 0, color = :black, ylims = [160, 260], yticks = [160, 260], size = [600, 150], shape = :x, xformatter = _ -> ""); png("temp")
# JLD2.@save "G:/BF/lorenz63/bfcnA.jld2" bfcn

@time _trajA0 = factory_lorenz63(DataFrame, [σ_[ 801], ρ_[ 801], b_[ 801]], saveat = 900:1e-4:1000)
@time _trajA1 = factory_lorenz63(DataFrame, [σ_[1201], ρ_[1201], b_[1201]], saveat = 900:1e-4:1000)
vrbl = reverse(half(names(_trajA0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:2)

for (η1, η2) = [[-1, -1], [-2, -2], [-1, -2], [-2, -1]]
    saveat = 1900:1e-3:2000
    @time trajA0 = add_diff(add_noise(_trajA0[:, last(vrbl)], exp10(η1)), method = :FDM, dt = 1e-5)
    @time trajA1 = add_diff(add_noise(_trajA1[:, last(vrbl)], exp10(η2)), method = :FDM, dt = 1e-5)

    f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f0 |> println
    f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f1 |> println
    trajB0 = ssolve(f0, trajA0[[1], f0.rname], saveat)
    trajB1 = ssolve(f1, trajA1[[1], f1.rname], saveat)

    # plot(
    #     plot(trajA0.x, trajA0.y, trajA0.z, alpha = .5, color = :black),
    #     plot(trajA1.x, trajA1.y, trajA1.z, alpha = .5, color = :black),
    #     plot(trajB0.x, trajB0.y, trajB0.z, alpha = .5, color = :red),
    #     plot(trajB1.x, trajB1.y, trajB1.z, alpha = .5, color = :red),
    #     ticks = [], layout = (:, 2), size = (400, 600)
    # )

    βm = (pm - p0) / (p1 - p0)
    βM = (pM - p0) / (p1 - p0)
    β0, β1 = 0, 1
    β_ = range(βm, βM, length = 2001)
    f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
    bfcn = callbfcn()
    @showprogress for k in eachindex(β_)
        sol = ssolve(f_[k], trajA0[[1], f_[k].rname], saveat)
        bfcn[β_[k]] = sol.z[arglmax(sol.z)]
    end
    # bfcn = callbfcn("G:/BF/lorenz63/bfcnB.jld2")
    # scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .1, ma = 1.0, msw = 0, color = :red, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> ""); png("temp")
    _bfcn = deepcopy(bfcn)
    bfcn = deepcopy(_bfcn)
    JLD2.@save "G:/BF/lorenz63/bfcnB_$(η1)$(η2).jld2" bfcn
end

for (η1, η2) = [[-1, -1], [-2, -2], [-1, -2], [-2, -1]]
    bfcn = JLD2.load("G:/BF/lorenz63/bfcnB_$(η1)$(η2).jld2")["bfcn"]
    scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .1, ma = 1.0, msw = 0, color = :red, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> ""); png("temp_$(η1)$(η2)")
end