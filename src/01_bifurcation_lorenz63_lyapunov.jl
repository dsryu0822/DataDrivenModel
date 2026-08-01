include.("../core/" .* readdir("core")[[1,2,3,4,6]])
using ChaosTools

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                Lorenz-63 ground truth

''''''''''''''''''''''''''''''''''''''''''''''''''"""
# function sys(du, u, p, t)
#     x, y, z = u; σ, ρ, β = p
    
#     du[1] = σ*(y - x)
#     du[2] = x*(ρ - z) - y
#     du[3] = x*y - β*z
#     return du
# end

# sol = factory_lorenz63(DataFrame, [10, 28, 8/3])
# plot(sol.x, sol.y, sol.z, alpha = .5)

pm, pM = -2, 3; p0, p1 = 0, 1;
p_ = range(pm, pM, length = 2001)
σ_ = range(6, 15, length = 2001)
ρ_ = range(120, 150, length = 2001)
b_ = range(3, 5, length = 2001)
# lpnv = callbfcn()
# @showprogress @threads for k in eachindex(p_)
#     lpnv[p_[k]] = lyapunovspectrum(CoupledODEs(sys, [100.0, 100, 100], [σ_[k], ρ_[k], b_[k]]), 100000)
# end
# JLD2.@save "G:/BF/lorenz63/lpnvA.jld2" lpnv

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                Lorenz-63 recovered

''''''''''''''''''''''''''''''''''''''''''''''''''"""
trajA0 = factory_lorenz63(DataFrame, [σ_[ 801], ρ_[ 801], b_[ 801]], saveat = 900:1e-3:1000)
trajA1 = factory_lorenz63(DataFrame, [σ_[1201], ρ_[1201], b_[1201]], saveat = 900:1e-3:1000)
vrbl = reverse(half(names(trajA0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:2)
f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f0 |> println
f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f1 |> println
trajB0 = ssolve(f0, trajA0[[1], f0.rname], 900:1e-3:1000)
trajB1 = ssolve(f1, trajA1[[1], f1.rname], 900:1e-3:1000)

# plot(
#     plot(trajA0.x, trajA0.y, trajA0.z, alpha = .5, color = :black),
#     plot(trajA1.x, trajA1.y, trajA1.z, alpha = .5, color = :black),
#     plot(trajB0.x, trajB0.y, trajB0.z, alpha = .5, color = :red),
#     plot(trajB1.x, trajB1.y, trajB1.z, alpha = .5, color = :red),
# )

βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 2001)
f_ = affine(Function, f0, f1) # affine(String, f0, f1) |> println
lpnv = callbfcn()
@showprogress @threads for k in eachindex(β_)
# for k in eachindex(β_)
    lpnv[β_[k]] = lyapunovspectrum(CoupledODEs(f_, [100.0, 100, 100], (β_[k],)), 10000)
end
# JLD2.@save "G:/BF/lorenz63/lpnvB.jld2" lpnv


# f_ = [define(Function, syntheticSINDy((1-β_[k])*f0.matrix + β_[k]*f1.matrix, vrbl, cnfg, method = "SINDy"), fname = "f_$(k)") for k in eachindex(β_)]
# lpnv = callbfcn()
# lpnv = JLD2.load("G:/BF/lorenz63/lpnvB.jld2")["lpnv"]
# for k in eachindex(β_)[401:end]
#     @time "k = $k" lpnv[β_[k]] = lyapunovspectrum(CoupledODEs(f_[k], [100.0, 100, 100], ()), 10000)
#     println("Live heap: ", Base.gc_live_bytes() / 1e6, " MB")
#     GC.gc()
#     if iszero(mod(k, 50))
#         @info "Saving lpnv to file at k = $k"
#         JLD2.@save "G:/BF/lorenz63/lpnvB.jld2" lpnv
#     end
# end
# JLD2.@save "G:/BF/lorenz63/lpnvB.jld2" lpnv

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                Lorenz-63 analyze

''''''''''''''''''''''''''''''''''''''''''''''''''"""
lpnvA = JLD2.load("G:/BF/lorenz63/lpnvA.jld2")["lpnv"]
df_lpnvA = sort(DataFrame([[keys(lpnvA)...] stack(values(lpnvA), dims = 1)], [:p, :λ1, :λ2, :λ3]), :p)
plt_lpnv = plot()
plot!(df_lpnvA.p, df_lpnvA.λ1, color = :black, label = "λ1")
plot!(df_lpnvA.p, df_lpnvA.λ2, color = :black, label = "λ2")
plot!(df_lpnvA.p, df_lpnvA.λ3, color = :black, label = "λ3")

lpnvB = JLD2.load("G:/BF/lorenz63/lpnvB.jld2")["lpnv"]
df_lpnvB = sort(DataFrame([[keys(lpnvB)...] stack(values(lpnvB), dims = 1)], [:β, :λ1, :λ2, :λ3]), :β)
plot!(df_lpnvB.β, df_lpnvB.λ1, color = :red, label = "λ1")
plot!(df_lpnvB.β, df_lpnvB.λ2, color = :red, label = "λ2")
plot!(df_lpnvB.β, df_lpnvB.λ3, color = :red, label = "λ3")
png("G:/lorenz_lyapunov.png")