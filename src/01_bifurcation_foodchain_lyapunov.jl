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


ic = [0.820915, 0.158239, 0.953786]
sol = solve(ODEProblem(f_, ic, (0, 10000)), p = [0])

βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 2001)
f_ = affine(Function, f0, f1) # affine(String, f0, f1) |> println
g_ = affine(Function, g0, g1) # affine(String, g0, g1) |> println
lpnv = callbfcn()
# for k in eachindex(β_)
@showprogress @threads for k in eachindex(β_)
    for l in 1:100        
        ic = collect(trajA0[1+100(l-1), 2:4])
        sol = solve(ODEProblem(f_, ic, (0, 10000)), p = [β_[k]], saveat = 0:1e-1:10000)
        if sol.retcode == ReturnCode.Success
            P_ = sol[3, sol.t .≥ 9000]
            if !isempty(P_) && minimum(P_) > 0.10
                lpnv[β_[k]] = lyapunovspectrum(CoupledODEs(f_, ic, (β_[k],)), 10000)
                break
            end
        end
    end
end
JLD2.@save "G:/BF/foodchain/lpnvB.jld2" lpnv

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                foodchain analyze

''''''''''''''''''''''''''''''''''''''''''''''"""
lpnvA = JLD2.load("G:/BF/foodchain/lpnvA.jld2")["lpnv"]
df_lpnvA = sort(DataFrame([[keys(lpnvA)...] stack(values(lpnvA), dims = 1)], [:p, :λ1, :λ2, :λ3]), :p)
lpnvB = JLD2.load("G:/BF/foodchain/lpnvB.jld2")["lpnv"]
df_lpnvB = sort(DataFrame([[keys(lpnvB)...] stack(values(lpnvB), dims = 1)], [:p, :λ1, :λ2, :λ3]), :p)
lpnvC = JLD2.load("G:/BF/foodchain/lpnvC.jld2")["lpnv"]
df_lpnvC = sort(DataFrame([[keys(lpnvC)...] stack(values(lpnvC), dims = 1)], [:β, :λ1, :λ2, :λ3]), :β)

plt_lpnv = plot(xlims = (pm, pM), xticks = [pm, p0, p1, pM])
plot!(df_lpnvA.p, df_lpnvA.λ1, color = :black, label = "λ1")
plot!(df_lpnvA.p, df_lpnvA.λ2, color = :black, label = "λ2")
plot!(df_lpnvA.p, df_lpnvA.λ3, color = :black, label = "λ3")
plot!(twiny(), df_lpnvB.p, df_lpnvB.λ1, color = :red, label = "λ1", xlims = (βm, βM), xticks = [0, 1])
plot!(twiny(), df_lpnvB.p, df_lpnvB.λ2, color = :red, label = "λ2", xlims = (βm, βM), xticks = [0, 1])
plot!(twiny(), df_lpnvB.p, df_lpnvB.λ3, color = :red, label = "λ3", xlims = (βm, βM), xticks = [0, 1])
plot!(twiny(), df_lpnvC.β, df_lpnvC.λ1, color = :blue, label = "λ1", xlims = (βm, βM), xticks = [0, 1])
plot!(twiny(), df_lpnvC.β, df_lpnvC.λ2, color = :blue, label = "λ2", xlims = (βm, βM), xticks = [0, 1])
plot!(twiny(), df_lpnvC.β, df_lpnvC.λ3, color = :blue, label = "λ3", xlims = (βm, βM), xticks = [0, 1])
png("G:/foodchain_lyapunov_sindyonly.png")


foo = callbfcn("G:/BF/foodchain/bfcnB_9596.jld2")

collect(keys(foo))[.!isempty.(values(foo))] |> sort


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    only four points

''''''''''''''''''''''''''''''''''''''''''''''''''"""


K = 0.96; k = 4
ic = [0.820915, 0.158239, 0.953786]
K_ = [0.93, 0.94, 0.95, 0.96]
f_ = []
g_ = []
ic_ = []
for k in eachindex(K_)
    K = K_[k]
    trajA = factory_foodchain(DataFrame, K, ic = ic, saveat = 4000:1e-1:5000)
    ic_ = push!(ic_, collect(trajA[1, 2:4]))

    vrbl = reverse(half(names(trajA[:, Not(:t)])))
    cnfg = cook(vrbl, poly = 0:4)
    f = define(Function, SINDy(trajA, vrbl, cnfg; λ = 1e-8), fname = "f$k");
    g = define(Function, SINDyPI(trajA, vrbl, cookPI(vrbl, poly = 0:3); λ = 1e-8), fname = "g$(rand(UInt8))");
    push!(f_, f)
    push!(g_, g)
end

results = []
for k in eachindex(K_)
    K = K_[k]
    println("K = $K")
    λo = lyapunovspectrum(CoupledODEs(sys, ic_[1], [K]), 1000000, Δt = 1e-1)
    println("λo = $λo")
    λf = lyapunovspectrum(CoupledODEs(f_[k], ic_[1], (1,)), 1000000, Δt = 1e-1)
    println("λf = $λf")
    λg = lyapunovspectrum(CoupledODEs(g_[k], ic_[1], (1,)), 1000000, Δt = 1e-1)
    println("λg = $λg")
    push!(results, [K; λo; λf; λg])
end

result = DataFrame((stack(results)')[:, [1,2,5,8,3,6,9,4,7,10]], [:K, :λo1, :λf1, :λg1, :λo2, :λf2, :λg2, :λo3, :λf3, :λg3])
CSV.write("lyapunov93949596.csv", result)

plot(
    bar(collect(result[3, 2:4]), color = [:black, :orangered, :dodgerblue]),
    bar(collect(result[3, 5:7]), color = [:black, :orangered, :dodgerblue]),
    bar(collect(result[3, 8:10]), color = [:black, :orangered, :dodgerblue]),
    layout = (1, 3)
)

define(String, f0, fname = "asdf") |> println


SINDy(trajA, vrbl, cnfg; λ = 1e-8) |> define |> println


λf = lyapunovspectrum(CoupledODEs(foo1, ic, (0,)), 20)

function foo1(dxyz, xyz, param, tau)
R, C, P = xyz; dR, dC, dP = dxyz;
dxyz[1] = -4.3418363406891345 + 10.657170746108159R + 9.809600036301244C + 7.887623691023628P - 7.024947615186845R*R - 17.77521055939072R*C - 16.455962882627635R*P - 7.010993907186865C*C - 14.849779217720897C*P - 3.963940692078413P*P + 0.6710754040109838R*R*R + 7.947888139968765R*R*C + 9.077727350040442R*R*P + 8.373889005955519R*C*C + 19.947901431837376R*C*P + 7.430931459873615R*P*P + 1.2827497416952498C*C*C + 7.961405931456596C*C*P + 5.572563190456326C*P*P + 0.3483261134703018P*P*P - 0.011388318836012018R*R*R*R - 0.8056964060080317R*R*R*C - 0.509958088747304R*R*R*P - 2.6322606649343387R*R*C*C - 5.210927801992969R*R*C*P - 3.445492892602658R*R*P*P - 0.9793755843787334R*C*C*C - 4.697067050852702R*C*C*P - 4.946824065883238R*C*P*P - 0.3652465553600266R*P*P*P - 0.07981006739528451C*C*C*C - 1.387222835861321C*C*C*P - 1.8145912617062203C*C*P*P - 0.3717502817837676C*P*P*P + 0.004880665792457781P*P*P*P
dxyz[2] = 3.7157373703298147 - 8.626297563433088R - 6.211696471289118C - 7.172564075479201P + 5.602243515410343R*R + 13.290483544561187R*C + 14.8862892353903R*P + 3.4229556754140558C*C + 7.859860972855601C*P + 4.419533274063009P*P - 0.6964815243928237R*R*R - 7.887162111288069R*R*C - 8.127725584103999R*R*P - 6.963425914532077R*C*C - 11.66208459038519R*C*P - 7.590396223018332R*P*P - 1.5662812354928692C*C*C - 1.7010117537643055C*C*P - 3.7172268390428056C*P*P - 0.882878038564709P*P*P + 0.013495819132554753R*R*R*R + 1.1528943720192197R*R*R*C + 0.45527206946731885R*R*R*P + 3.014344556862239R*R*C*C + 4.218309614182419R*R*C*P + 3.0402339556529183R*R*P*P + 1.343110750327052R*C*C*C + 2.8208414328334324R*C*C*P + 1.7487677568561433R*C*P*P + 0.9966891213984526R*P*P*P + 0.23916903841386433C*C*C*C + 1.4454742302515196C*C*C*P - 0.5940847060168469C*C*P*P + 1.1546877994511178C*P*P*P - 0.03981454661465849P*P*P*P
dxyz[3] = 0.6260989686557448 - 1.0308731799045094R - 3.9979035597076233C - 0.7950596121049894P + 0.38103743217346303R*R + 4.4847270095132625R*C + 1.5696736413907866R*P + 3.5880382287757984C*C + 6.989918234870334C*P - 0.4555925836163925P*P + 0.02540612027434092R*R*R - 0.06072602879001906R*R*C - 0.9500017636299101R*R*P - 1.410463090825943R*C*C - 8.28581683078979R*C*P + 0.1594647660240564R*P*P + 0.2835314929704764C*C*C - 6.26039417224586C*C*P - 1.8553363470399913C*P*P + 0.5345519249735269P*P*P - 0.0021075002917422467R*R*R*R - 0.3471979657427143R*R*R*C + 0.05468601932552348R*R*R*P - 0.3820838917993906R*R*C*C + 0.9926181870749483R*R*C*P + 0.40525893572902066R*R*P*P - 0.36373516553598256R*C*C*C + 1.8762256170517988R*C*C*P + 3.1980563043217396R*C*P*P - 0.6314425659134127R*P*P*P - 0.159358971112489C*C*C*C - 0.058251393412941777C*C*C*P + 2.408675965020895C*C*P*P - 0.7829375173133204C*P*P*P + 0.034933880816479704P*P*P*P
end

sol = solve(ODEProblem(foo1, ic, (0, 10000)))