include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    food chain

''''''''''''''''''''''''''''''''''''''''''''''''''"""
# sol = factory_foodchain(DataFrame, .99, ic = [.85, .2, .8], saveat = 0:1e-1:10000)
# plot(sol.R, sol.C, sol.P)

pm, pM = .88, 1.00; p0, p1 = .940, .950;
K_ = range(pm, pM, length = 1001)
# bfcn = callbfcn()
# @showprogress @threads for k in eachindex(K_)
#     for _ in 1:100
#         sol = factory_foodchain(DataFrame, K_[k], ic = [.85, 0.5rand(), .8], saveat = 0:1e-1:10000)
#         P_ = sol.P[sol.t .≥ 9000]
#         if !isempty(P_) && minimum(P_) > 0.55
#             bfcn[K_[k]] = P_[arglmin(P_)]
#             break
#         end
#     end
# end
# # bfcn = callbfcn("G:/BF/foodchain/bfcnA.jld2")
# scatter(dict2bifurcation(bfcn)..., xticks = [pm, p0, p1, pM], ms = .5, ma = .5, msw = 0, color = :black); png("temp")
# # JLD2.@save "G:/BF/foodchain/bfcnA.jld2" bfcn

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

# plot(
#     plot(trajA0.R, trajA0.C, trajA0.P, alpha = .5, color = :black),
#     plot(trajA1.R, trajA1.C, trajA1.P, alpha = .5, color = :black),
#     plot(trajB0.R, trajB0.C, trajB0.P, alpha = .5, color = :red),
#     plot(trajB1.R, trajB1.C, trajB1.P, alpha = .5, color = :red),
#     plot(trajC0.R, trajC0.C, trajC0.P, alpha = .5, color = :blue),
#     plot(trajC1.R, trajC1.C, trajC1.P, alpha = .5, color = :blue),
#     layout = (:, 2), size = [400, 600]
# )
# CSV.write("G:/BF/foodchain/trajA0_9394.csv", trajA0)
# CSV.write("G:/BF/foodchain/trajA1_9394.csv", trajA1)
# CSV.write("G:/BF/foodchain/trajB0_9394.csv", trajB0)
# CSV.write("G:/BF/foodchain/trajB1_9394.csv", trajB1)
# CSV.write("G:/BF/foodchain/trajC0_9394.csv", trajC0)
# CSV.write("G:/BF/foodchain/trajC1_9394.csv", trajC1)

βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 1001)
cnfg = cook(vrbl, poly = 0:4)
f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
bfcn = callbfcn()
@showprogress @threads for k in eachindex(β_)
    for _ in 1:10
        sol = ssolve(f_[k], trajA0[[1], f_[k].rname], 0:1e-1:10000)
        P_ = sol.P[sol.t .≥ 9000]
        if !isempty(P_) && minimum(P_) > 0.55
            bfcn[β_[k]] = P_[arglmin(P_)]
            break
        end
    end
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], ms = .5, ma = .5, msw = 0, color = :red, size = [400, 200], xlims = [βm, βM]); png("temp1")
JLD2.@save "G:/BF/foodchain/bfcnB_9495.jld2" bfcn

cnfg = cookPI(vrbl, poly = 0:3)
g_ = [syntheticSINDy((1-β)*g0.matrix + β*g1.matrix, vrbl, cnfg, method = "SINDyPI") for β in β_]
bfcn = callbfcn("G:/BF/foodchain/bfcnC.jld2")
@showprogress @threads for k in eachindex(β_)
    for _ in 1:10
        sol = ssolve(g_[k], trajA0[[1], g_[k].rname], 0:1e-1:10000)
        P_ = sol.P[sol.t .≥ 9000]
        if !isempty(P_) && minimum(P_) > 0.55
            bfcn[β_[k]] = P_[arglmin(P_)]
            break
        end
    end
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .5, msw = 0, color = :blue); png("temp2")
JLD2.@save "G:/BF/foodchain/bfcnC_9495.jld2" bfcn

# """''''''''''''''''''''''''''''''''''''''''''''''''''

#             pitchfork bifurcation

# ''''''''''''''''''''''''''''''''''''''''''''''''''"""
# sol = factory_foodchain(DataFrame, 0.9612, ic = [.85, .4, .8], saveat = 9000:1e-1:9500)
# ic1 = collect(sol[argmin(sol.P), [:R, :C, :P]])
# sol = factory_foodchain(DataFrame, 0.9612, ic = [.85, .2, .8], saveat = 9000:1e-1:9500)
# ic2 = collect(sol[argmin(sol.P), [:R, :C, :P]])

# # K_ = 0.96080:1e-5:0.96150
# K_ = [0.96, 0.960802, 0.961, 0.961531, 0.962]
# plt_ = []
# @showprogress for k in eachindex(K_)
#     plt = plot(legend = :none, framestyle = :none, ticks = [])
#     sol = factory_foodchain(DataFrame, K_[k], ic = ic1, saveat = 0:1e-1:3000)[(end-1000):end, :]
#     CSV.write("G:/BF/foodchain/trajA_1$(k).csv", sol)
#     plot!(plt, sol.R, sol.C, sol.P; color = 1, lw = 2, alpha = .5)
#     sol = factory_foodchain(DataFrame, K_[k], ic = ic2, saveat = 0:1e-1:3000)[(end-1000):end, :]
#     CSV.write("G:/BF/foodchain/trajA_2$(k).csv", sol)
#     plot!(plt, sol.R, sol.C, sol.P; color = 2, lw = 2, alpha = .5)
#     push!(plt_, plt)
#     # png(plt, "G:/K_$(rpad(K_[k], 8, '0')).png")
# end
# plot(plt_..., layout = (1, :), size = (800, 200), margin = - 0mm)

# pm, pM = .9606, 0.9616; p0, p1 = .950, .960;
# K_ = range(pm, pM, length = 1001)
# bfcn1 = callbfcn()
# bfcn2 = callbfcn()
# @showprogress @threads for k in eachindex(K_)
#     sol1 = factory_foodchain(DataFrame, K_[k], ic = ic1, saveat = 0:1e-1:10000)
#     sol2 = factory_foodchain(DataFrame, K_[k], ic = ic2, saveat = 0:1e-1:10000)
#     z1_ = sol1.P[sol1.t .≥ 9000]
#     z2_ = sol2.P[sol2.t .≥ 9000]
#     bfcn1[K_[k]] = z1_[arglmin(z1_)]
#     bfcn2[K_[k]] = z2_[arglmin(z2_)]
# end
# # JLD2.@save "G:/BF/foodchain/bfcnA_1_9606_9616.jld2" bfcn1
# # JLD2.@save "G:/BF/foodchain/bfcnA_2_9606_9616.jld2" bfcn2
# bfcn1 = callbfcn("G:/BF/foodchain/bfcnA_1_9606_9616.jld2")
# bfcn2 = callbfcn("G:/BF/foodchain/bfcnA_2_9606_9616.jld2")
# plot()
# scatter!(dict2bifurcation(bfcn1)..., xticks = [pm, p0, p1, pM], ms = .5, ma = .5, msw = 0, color = 1);
# scatter!(dict2bifurcation(bfcn2)..., xticks = [pm, p0, p1, pM], ms = .5, ma = .5, msw = 0, color = 2); png("temp")



# """''''''''''''''''''''''''''''''''''''''''''''''''''

#             different initial conditions

# ''''''''''''''''''''''''''''''''''''''''''''''''''"""

# # ic_ = [collect(sol[k, [:R, :C, :P]]) for k in shuffle(1:nrow(sol))[1:100]]
# ic_ = DataFrame(R = rand(100), C = rand(100), P = rand(100))
# plt_attA = plot(legend = :none)
# plt_attB = plot(legend = :none)
# plt_attC = plot(legend = :none)
# @showprogress for (k, ic) = enumerate(eachrow(ic_))
#     trajA = factory_foodchain(DataFrame, p1, ic = [ic...], saveat = 900:1e-1:1000)
#     trajB = ssolve(f1, ic, 900:1e-1:1000)
#     trajC = ssolve(g1, ic, 900:1e-1:1000)
#     # CSV.write("G:/BF/foodchain/tempA95_$k.csv", trajA)
#     # CSV.write("G:/BF/foodchain/tempB95_$k.csv", trajB)
#     # CSV.write("G:/BF/foodchain/tempC95_$k.csv", trajC)
#     plot!(plt_attA, trajA.R, trajA.C, trajA.P; color = :black)
#     plot!(plt_attB, trajB.R, trajB.C, trajB.P; color = :red)
#     plot!(plt_attC, trajC.R, trajC.C, trajC.P; color = :blue)
# end
# plot(plt_attA, plt_attB, plt_attC; layout = (1, 3), size = (1200, 400)); png("temp")

# ic_ = DataFrame(stack([Base.product([repeat(0:1e-1:1) for _ in 1:3]...)...], dims = 1), [:R, :C, :P])
# CSV.write("G:/BF/foodchain/vtfdA95.csv", DataFrame(first.([factory_foodchain(DataFrame, p0, ic = [ic...], saveat = 0:1e-1:1) for ic in eachrow(ic_)]))[:, Not(:t)])
# CSV.write("G:/BF/foodchain/vtfdB95.csv", [ic_ f0(ic_)])
# CSV.write("G:/BF/foodchain/vtfdC95.csv", [ic_ g0(ic_)])