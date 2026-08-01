include.("../core/" .* readdir("core")[[1,2,3,4,6]])
using ChaosTools

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                food chain ground truth

''''''''''''''''''''''''''''''''''''''''''''''''''"""
function sys(du, u, p, t)
    xc,    yc,   xp,    yp,      R0,  C0 = (
0.4, 2.009, 0.08, 2.876, 0.16129, 0.5)
    R,C,P = u
    K = p[1]

    du[1] = R*(1 - frac(R,K)) - xc*yc*frac(C*R, R + R0)
    du[2] = xc*C*(frac(yc*R, R + R0) - 1) - xp*yp*frac(P*C, C + C0)
    du[3] = xp*P*(frac(yp*C, C + C0) - 1)
    return du
end
# sol = factory_foodchain(DataFrame, .99, ic = [.85, .2, .8], saveat = 0:1e-1:10000)
# plot(sol.R, sol.C, sol.P)


pm, pM = .88, 1.00; p0, p1 = 0.95, 0.96;
K_ = range(pm, pM, length = 2001)
lpnv = callbfcn()
@showprogress @threads for k in eachindex(K_)
    for _ in 1:10
        ic = [.85, 0.5rand(), .8]
        sol = factory_foodchain(DataFrame, K_[k], ic = ic, saveat = 0:1e-1:10000)
        P_ = sol.P[sol.t .≥ 9000]
        if !isempty(P_) && minimum(P_) > 0.55
            lpnv[K_[k]] = lyapunovspectrum(CoupledODEs(sys, ic, [K_[k]]), 10000)
            break
        end
    end
end
JLD2.@save "G:/BF/foodchain/lpnvA.jld2" lpnv

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                food chain recovered

''''''''''''''''''''''''''''''''''''''''''''''''''"""
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

βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 2001)
g_ = affine(Function, g0, g1) # affine(String, g0, g1) |> println
lpnv = callbfcn()
@showprogress @threads for k in eachindex(β_)
    for _ in 1:10
        ic = [.85, 0.5rand(), .8]
        sol = ssolve(f_, ic, 4000:1e-2:5000, (β_[k],))
        P_ = sol.P[sol.t .≥ 4900]
        if !isempty(P_) && minimum(P_) > 0.55
            lpnv[β_[k]] = lyapunovspectrum(CoupledODEs(f_, ic, (β_[k],)), 10000)
            break
        end
    end
end
JLD2.@save "G:/BF/foodchain/lpnvB.jld2" lpnv