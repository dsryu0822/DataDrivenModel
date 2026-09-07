include.("../core/" .* readdir("core")[[1,2,3,4,6]])

ReLU(x) = max(0, x)
function factory_foodchain_noisy(K::Number, ξ = 0; ic = [.85, .12 + .34rand(), .8], saveat = 0:1e-2:10, seed = -1, kargs...)
     xc,    yc,   xp,    yp,      R0,  C0 = (
    0.4, 2.009, 0.08, 2.876, 0.16129, 0.5)
    dt = saveat.step.hi
    if seed > 0 Random.seed!(seed) end
    function sys(du, u, p, t)
        R,C,P = ReLU.(u)
        K = p[1]

        du[1] = R*(1 - frac(R,K)) - xc*yc*frac(C*R, R + R0) + ξ*sqrt(R*dt)*randn()
        du[2] = xc*C*(frac(yc*R, R + R0) - 1) - xp*yp*frac(P*C, C + C0) + ξ*sqrt(C*dt)*randn()
        du[3] = xp*P*(frac(yp*C, C + C0) - 1) + ξ*sqrt(P*dt)*randn()
        return du
    end
    sol = solve(ODEProblem(sys, ic, (0, last(saveat)), [K]), Heun(), dt = saveat.step.hi, adaptive=false, maxiters = Inf)
    dsol = [sol.k[2][1] stack([solk[2] for solk in sol.k[2:end]])]
    matrix = [sol.t sol[:, :]' dsol']
    return matrix[sol.t .≥ first(saveat), :]
end
factory_foodchain_noisy(T::Type, args...; kargs...) =
DataFrame(factory_foodchain_noisy(args...; kargs...), ["t", "R", "C", "P", "dR", "dC", "dP"])

# K_ = 0.93:1e-2:0.96
# ξ_ = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1e-0]
# plt_ = []
# for K in K_
#     for ξ in ξ_
#         Random.seed!(4) # seed = [1,2,6]
#         sol = factory_foodchain_noisy(DataFrame, K, ξ, ic = [0.820915, 0.158239, 0.953786], saveat = 9000:1e-1:10000)
#         push!(plt_, plot(sol.R, sol.C, sol.P,
#         lw = 0.5, alpha = 0.5,
#         color = :black, legend = :none, formatter = _ -> ""))
#     end
# end

# plot(plt_..., layout = (length(K_), :), size = (length(ξ_)*200, length(K_)*200), camera = [75, 30], dpi = 180); png("temp")

g93 = SINDyPI(trajA93, vrbl, cnfg; λ = 1e-8); g93 |> println
cnfg = (cnfg[1][.!iszero.(sum.(eachrow(g93.matrix))), :], cnfg[2][.!iszero.(sum.(eachrow(g93.matrix))), :])

ξ = 1e-2
λ = 1e-8
seed = 16
ic0 = [0.820915, 0.158239, 0.953786]
trajA93 = factory_foodchain_noisy(DataFrame, 0.93, ξ, ic = ic0, saveat = 9000:1e-1:10000, seed = seed); vrbl = reverse(half(names(trajA93[:, Not(:t)])))
trajA94 = factory_foodchain_noisy(DataFrame, 0.94, ξ, ic = ic0, saveat = 9000:1e-1:10000, seed = seed)
trajA95 = factory_foodchain_noisy(DataFrame, 0.95, ξ, ic = ic0, saveat = 9000:1e-1:10000, seed = seed)
trajA96 = factory_foodchain_noisy(DataFrame, 0.96, ξ, ic = ic0, saveat = 9000:1e-1:10000, seed = seed)

cnfg = cook(vrbl, poly = 0:4)
f93 = SINDy(trajA93, vrbl, cnfg; λ); f93 |> println
f94 = SINDy(trajA94, vrbl, cnfg; λ); f94 |> println
f95 = SINDy(trajA95, vrbl, cnfg; λ); f95 |> println
f96 = SINDy(trajA96, vrbl, cnfg; λ); f96 |> println

trajB93 = ssolve(f93, trajA93[[1], f93.rname], 4000:1e-1:5000)
trajB94 = ssolve(f94, trajA94[[1], f94.rname], 4000:1e-1:5000)
trajB95 = ssolve(f95, trajA95[[1], f95.rname], 4000:1e-1:5000)
trajB96 = ssolve(f96, trajA96[[1], f96.rname], 4000:1e-1:5000)

plot(
    plot(trajA93.R, trajA93.C, trajA93.P, alpha = .5, color = :black),
    plot(trajA94.R, trajA94.C, trajA94.P, alpha = .5, color = :black),
    plot(trajA95.R, trajA95.C, trajA95.P, alpha = .5, color = :black),
    plot(trajA96.R, trajA96.C, trajA96.P, alpha = .5, color = :black),
    plot(trajB93.R, trajB93.C, trajB93.P, alpha = .5, color = :red),
    plot(trajB94.R, trajB94.C, trajB94.P, alpha = .5, color = :red),
    plot(trajB95.R, trajB95.C, trajB95.P, alpha = .5, color = :red),
    plot(trajB96.R, trajB96.C, trajB96.P, alpha = .5, color = :red),
    layout = (2, :), size = [800, 400], formatter = _ -> ""
); png("temp")


cnfg = cookPI(vrbl, poly = 0:3)
# cnfg = (g93.recipeF, cnfg[2][g93.recipeF.index, :])
# cnfg = cookPI(vrbl, poly = 0:3)[[2,3,4,5,6,8,9,11,14,15,18,22,42,43,46,63,68], :]
g93 = SINDyPI(trajA93, vrbl, cnfg; λ); g93 |> println
g94 = SINDyPI(trajA94, vrbl, cnfg; λ); g94 |> println
g95 = SINDyPI(trajA95, vrbl, cnfg; λ); g95 |> println
g96 = SINDyPI(trajA96, vrbl, cnfg; λ); g96 |> println

trajC93 = ssolve(g93, trajA93[[1], g93.rname], 4000:1e-1:5000)
trajC94 = ssolve(g94, trajA94[[1], g94.rname], 4000:1e-1:5000)
trajC95 = ssolve(g95, trajA95[[1], g95.rname], 4000:1e-1:5000)
trajC96 = ssolve(g96, trajA96[[1], g96.rname], 4000:1e-1:5000)

plot(
    plot(trajA93.R, trajA93.C, trajA93.P, alpha = .5, color = :black),
    plot(trajA94.R, trajA94.C, trajA94.P, alpha = .5, color = :black),
    plot(trajA95.R, trajA95.C, trajA95.P, alpha = .5, color = :black),
    plot(trajA96.R, trajA96.C, trajA96.P, alpha = .5, color = :black),
    plot(trajC93.R, trajC93.C, trajC93.P, alpha = .5, color = :blue),
    plot(trajC94.R, trajC94.C, trajC94.P, alpha = .5, color = :blue),
    plot(trajC95.R, trajC95.C, trajC95.P, alpha = .5, color = :blue),
    plot(trajC96.R, trajC96.C, trajC96.P, alpha = .5, color = :blue),
    layout = (2, :), size = [800, 400], formatter = _ -> ""
); png("temp")



trajA931 = factory_foodchain_noisy(DataFrame, 0.93, 0, ic = [0.820915, 0.158239, 0.953786], saveat = 9000:1e-1:10000);
trajA932 = factory_foodchain_noisy(DataFrame, 0.93, 1e-3, ic = [0.820915, 0.158239, 0.953786], saveat = 9000:1e-1:10000);
rmse(Matrix(trajA931[:, Not(:t)]), Matrix(trajA932[:, Not(:t)]))

g931 = SINDyPI(trajA931, vrbl, cnfg; λ = 1e-8); g931 |> println
# cnfg = (cnfg[1][.!iszero.(sum.(eachrow(g931.matrix))), :], cnfg[2][.!iszero.(sum.(eachrow(g931.matrix))), :])
g932 = SINDyPI(trajA932, vrbl, cnfg; λ = 1e-4); g932 |> println
g932 = SINDyPI(trajA932, vrbl, cnfg; λ = 1e-2); g932 |> println



seed = 6
for ξ = [-5, -4, -3, -2]
    saveat = 4000:1e-1:5000
# for (P0, P1) = [[93, 94], [93, 95], [94, 95], [93, 96], [94, 96], [95, 96]]
for (P0, P1) = [[94, 95]]
    # P0 = 95; P1 = 96
    pm = 0.88; pM = 1.00
    p0, p1 = P0/100, P1/100
    @info "G:/BF/foodchain/$(ξ)/bfcnB_$(P0)$(P1).jld2"

    βm = (pm - p0) / (p1 - p0)
    βM = (pM - p0) / (p1 - p0)
    β0, β1 = 0, 1
    β_ = range(βm, βM, length = 2001)

    trajA0 = factory_foodchain_noisy(DataFrame, p0, exp10(ξ), ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000, seed = seed)
    trajA1 = factory_foodchain_noisy(DataFrame, p1, exp10(ξ), ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000, seed = seed)
    vrbl = reverse(half(names(trajA0[:, Not(:t)])))
    cnfg = cook(vrbl, poly = 0:4)
    f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-8); f0 |> println
    f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-8); f1 |> println

    f_ = affine(Function, f0, f1)
    bfcn = callbfcn()
    @showprogress @threads for k in eachindex(β_)
        for l in 1:10
            ic = collect(trajA1[1+100(l-1), 2:4])
            sol = solve(ODEProblem(f_, ic, (0, 5000)), p = [β_[k]], RK4(), dt = saveat.step.hi, adaptive=false, maxiters = Inf)
            if sol.retcode == ReturnCode.Success
                P_ = sol[3, sol.t .≥ 4000]
                if !isempty(P_) && (minimum(P_) > 0.10)
                    bfcn[β_[k]] = P_[arglmin(P_)]
                    break
                end
            end
        end
    end
    # scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], ms = .5, ma = .5, msw = 0, color = :red, size = [400, 200], xlims = [βm, βM]); png("temp1")
    _bfcn = deepcopy(bfcn)
    bfcn = deepcopy(_bfcn)
    JLD2.@save "G:/BF/foodchain/$(ξ)/bfcnB_$(P0)$(P1).jld2" bfcn
end
end

β0, β1 = 0, 1
pm, pM = 0.88, 1.0
bfcnAh, bfcnAv = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnA.jld2"));
plt_ = []
for ξ = string.([-5, -4, -3, -2])
    (P0, P1) = (94, 95)
    @info "G:/BF/foodchain/$(ξ)/bfcnB_$(P0)$(P1).jld2"
    bfcnh, bfcnv = dict2bifurcation(callbfcn("G:/BF/foodchain/$(ξ)/bfcnB_$(P0)$(P1).jld2"))
    p0 = P0/100; p1 = P1/100;
    βm = (pm - p0) / (p1 - p0)
    βM = (pM - p0) / (p1 - p0)
    plt = scatter(bfcnAh, bfcnAv, color = :black, msw = 0, ms = 0.3, ylims = [0.55, 0.8], yticks = ([0.55, 0.8], ["0.55", "0.8"]), xticks = ([pm, p0, p1, pM], tickfmt.([βm, β0, β1, βM])), xlims = [pm, pM]);
    annotate!(plt, (-0.1, 0.5), text(L"P", 12, rotation = 90))
    annotate!(plt, (0.25, -0.15), text(L"\beta", 12))
    annotate!(plt, (0.25, 1.10), text(L"K", 11))
    scatter!(twiny(plt), bfcnh, bfcnv, color = :red, msw = 0, ms = 0.3, ylims = [0.55, 0.8], xticks = ([β0, β1], ["$p0   ", "   $p1"]), xlims = [βm, βM], grid = true);
    scatter!(twinx(plt), yticks = []); 
    push!(plt_, plt)
end
plot(plt_..., layout = (:, 1), size = [600, 900]); @time png("9495")

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    bistability

''''''''''''''''''''''''''''''''''''''''''''''''''"""

ic_ = DataFrame(R = rand(20), C = rand(20), P = rand(20))
# for seed = 631:90
    plt_att94_ = []
    plt_att95_ = []
for ξ = [1e-5, 1e-4, 1e-3, 1e-2]
    plt_att94 = plot(legend = :none)
    plt_att95 = plot(legend = :none)
    trajA94 = factory_foodchain_noisy(DataFrame, 0.94, ξ, ic = ic0, saveat = 4000:1e-1:5000, seed = seed)
    trajA95 = factory_foodchain_noisy(DataFrame, 0.95, ξ, ic = ic0, saveat = 4000:1e-1:5000, seed = seed)
    vrbl = reverse(half(names(trajA94[:, Not(:t)])))
    cnfg = cook(vrbl, poly = 0:4)
    f94 = SINDy(trajA94, vrbl, cnfg; λ)
    f95 = SINDy(trajA95, vrbl, cnfg; λ)
    @showprogress for ic = eachrow(ic_)
        trajB94 = ssolve(f94, ic, 4000:1e-1:5000)
        trajB95 = ssolve(f95, ic, 4000:1e-1:5000)
        plot!(plt_att94, trajB94.R, trajB94.C, trajB94.P; color = 1, lims = [0, 1.2], formatter = _ -> "", ticks = [0, 1.2])
        plot!(plt_att95, trajB95.R, trajB95.C, trajB95.P; color = 2, lims = [0, 1.2], formatter = _ -> "", ticks = [0, 1.2])
    end
    push!(plt_att94_, plt_att94)
    push!(plt_att95_, plt_att95)
end
plot(
    plt_att94_[1], plt_att95_[1],
    plt_att94_[2], plt_att95_[2],
    plt_att94_[3], plt_att95_[3],
    plt_att94_[4], plt_att95_[4],
    lims = [0, 1.2], formatter = _ -> "",
    layout = (:, 2), size = 300 .* [2, 4], ticks = [0, 1.2],
    margin = -2mm
# ); @time "✅png" png("temp3")
); @time "✅png" png("temp3")
# end

lots.plotly()
Plots.gr()


plot(
    plt_[1], plt_att94_[1], plt_att95_[1],
    plt_[2], plt_att94_[2], plt_att95_[2],
    plt_[3], plt_att94_[3], plt_att95_[3],
    plt_[4], plt_att94_[4], plt_att95_[4],
    layout = grid(4, 3, widths = [0.6 ,0.2, 0.2]),
    size = 300 .* [3, 3],
    margin = -2mm
); @time "✅png" png("temp")
