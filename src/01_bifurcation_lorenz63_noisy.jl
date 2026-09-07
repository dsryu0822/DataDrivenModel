include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Lorenz-63

''''''''''''''''''''''''''''''''''''''''''''''''''"""
pm, pM = -2, 3; p0, p1 = 0, 1;
p_ = range(pm, pM, length = 2001)
σ_ = range(6, 15, length = 2001)
ρ_ = range(120, 150, length = 2001)
b_ = range(3, 5, length = 2001)
βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 2001)
# β_ = range(-100, 100, length = 2001)


@time _trajA0 = factory_lorenz63(DataFrame, [σ_[ 801], ρ_[ 801], b_[ 801]], saveat = 900:1e-3:1000)
@time _trajA1 = factory_lorenz63(DataFrame, [σ_[1201], ρ_[1201], b_[1201]], saveat = 900:1e-3:1000)
vrbl = reverse(half(names(_trajA0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:2)

for (η1, η2) = [[-3, -3], [-2, -2], [-1, -2], [-2, -1]] # η1, η2 = [-1, -1]
    # for (η1, η2) = [[-1, -1], [-0, -0], [-1, -0], [-0, -1]] # η1, η2 = [-0, -0]
    Random.seed!(1)
    saveat = 1900:1e-3:2000
    trajA0 = deepcopy(_trajA0[1:1:end, :]); trajA0.x .+= exp10(η1)*randn(nrow(trajA0)); trajA0.y .+= exp10(η1)*randn(nrow(trajA0)); trajA0.z .+= exp10(η1)*randn(nrow(trajA0))
    trajA1 = deepcopy(_trajA1[1:1:end, :]); trajA1.x .+= exp10(η2)*randn(nrow(trajA1)); trajA1.y .+= exp10(η2)*randn(nrow(trajA1)); trajA1.z .+= exp10(η2)*randn(nrow(trajA1))
    # trajA0 = shuffle(trajA0)[1:1000, :]; trajA1 = shuffle(trajA1)[1:1000, :]
    # plt_data0 = plot(_trajA0.x, _trajA0.y, _trajA0.z, alpha = .5, color = :black)
    # scatter!(plt_data0, trajA0.x, trajA0.y, trajA0.z, color = 2, msw = 0, ms = 2)
    # plt_data1 = plot(_trajA1.x, _trajA1.y, _trajA1.z, alpha = .5, color = :black)
    # scatter!(plt_data1, trajA1.x, trajA1.y, trajA1.z, color = 2, msw = 0, ms = 2)
    # plot(plt_data0, plt_data1, ticks = [], layout = (:, 2), size = (400, 200)); png("temp")

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

    # β_ = range(βm, βM, length = 2001)
    # f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
    f_ = affine(Function, f0, f1) # affine(String, f0, f1) |> println
    bfcn = callbfcn()
    @showprogress @threads for k in eachindex(β_)
        sol = ssolve(f_[k], trajA0[[1], f_[k].rname], saveat)
        bfcn[β_[k]] = sol.z[arglmax(sol.z)]
    end
    # bfcn = callbfcn("G:/BF/lorenz63/bfcnB.jld2")
    # scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .1, ma = 1.0, msw = 0, color = :red, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> ""); png("temp")
    _bfcn = deepcopy(bfcn)
    bfcn = deepcopy(_bfcn)
    JLD2.@save "G:/BF/lorenz63/bfcnB_$(η1)$(η2).jld2" bfcn
end


bfcnA = callbfcn("G:/BF/lorenz63/bfcnA.jld2")
plt_original = scatter(dict2bifurcation(bfcnA)..., xlims = [pm, pM], xticks = [pm, p0, p1, pM], ms = .1, ma = 1.0, msw = 0, color = :black, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> "");
# annotate!(plt_original, (0.05, 1.1), "(a)");
# annotate!(plt_original, (.5, -0.1), L"\mu");
# annotate!(plt_original, (-0.1, .5), text(L"z", rotation = 90));


plt_bifcn = [plt_original, plt_η]

η_ = [[-3, -3], [-2, -1], [-2, -2], [-1, -2], [-1, -1], [-1, 0], [0, 0], [0, -1]]
for (k, (η1, η2)) = enumerate(η_)
    @info "η1 = $η1, η2 = $η2"
    bfcn = JLD2.load("G:/BF/lorenz63/bfcnB_$(η1)$(η2).jld2")["bfcn"]
    plt_cases = scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .1, ma = 1.0, msw = 0, color = :red, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> "")
    # annotate!(plt_cases, (0.02, 1.1), "$(('c':'j')[k])");
    # annotate!(plt_cases, (.5, -0.1), L"\beta");
    # annotate!(plt_cases, (-0.04, .5), text(L"z", rotation = 90)); png("temp")
    push!(plt_bifcn, plt_cases)
end
plot(
    plt_bifcn[1], plt_bifcn[3], plt_bifcn[5], plt_bifcn[7], plt_bifcn[9], 
    plt_bifcn[2], plt_bifcn[4], plt_bifcn[6], plt_bifcn[8], plt_bifcn[10]
    , layout = (2, :), size = [1200, 400], bottom_margin = 2mm, top_margin = 10mm
); @time png("temp")




bfcnA = callbfcn("G:/BF/lorenz63/bfcnA.jld2")
# plt_original = scatter(dict2bifurcation(bfcnA)..., xlims = [pm, pM], xticks = [pm, p0, p1, pM], ms = .1, ma = 1.0, msw = 0, color = :black, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> "");
η_ = [[-2, -2], [-1, -1], [0, 0]]
for (k, (η1, η2)) = enumerate(η_)
    @info "η1 = $η1, η2 = $η2"
    bfcn = JLD2.load("G:/BF/lorenz63/bfcnB_$(η1)$(η2).jld2")["bfcn"]
    plt_original = scatter(dict2bifurcation(bfcnA)..., xlims = [pm, pM], xticks = [pm, p0, p1, pM], ms = .1, ma = 1.0, msw = 0, color = :darkgray, ylims = [160, 260], yticks = [160, 260], size = [600, 150], xformatter = _ -> "");
    plt_cases = scatter!(twiny(plt_original), dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .1, ma = 1.0, msw = 0, color = :blue, ylims = [160, 260], yticks = [160, 260], size = [600, 150]);
    @time png("temp$k")
end

qwer = plot(rand(10))
asdf = scatter!(twiny(qwer), 100:109, rand(10))


function lorenz_matrix(σ, ρ, b)
    A = sparse(zeros(10, 3))
    A[2, 1] = -σ
    A[3, 1] = σ
    A[2, 2] = ρ
    A[3, 2] = -1
    A[7, 2] = -1
    A[4, 3] = -b
    A[6, 3] = 1
    return A
end


trajA0 = deepcopy(_trajA0)
trajA1 = deepcopy(_trajA1)
f01 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f01 |> println
f11 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f11 |> println
define(f01) |> println
define(f11) |> println

affine(String, f01, f11) |> println

η1, η2 = [-1, -1]
Random.seed!(1)
saveat = 1900:1e-3:2000
trajA0 = deepcopy(_trajA0[1:100:end, :]); trajA0.x .+= exp10(η1)*randn(nrow(trajA0)); trajA0.y .+= exp10(η1)*randn(nrow(trajA0)); trajA0.z .+= exp10(η1)*randn(nrow(trajA0))
trajA1 = deepcopy(_trajA1[1:100:end, :]); trajA1.x .+= exp10(η2)*randn(nrow(trajA1)); trajA1.y .+= exp10(η2)*randn(nrow(trajA1)); trajA1.z .+= exp10(η2)*randn(nrow(trajA1))
f02 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f02 |> println
f12 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f12 |> println

η1, η2 = [0, 0]
Random.seed!(1)
saveat = 1900:1e-3:2000
trajA0 = deepcopy(_trajA0[1:100:end, :]); trajA0.x .+= exp10(η1)*randn(nrow(trajA0)); trajA0.y .+= exp10(η1)*randn(nrow(trajA0)); trajA0.z .+= exp10(η1)*randn(nrow(trajA0))
trajA1 = deepcopy(_trajA1[1:100:end, :]); trajA1.x .+= exp10(η2)*randn(nrow(trajA1)); trajA1.y .+= exp10(η2)*randn(nrow(trajA1)); trajA1.z .+= exp10(η2)*randn(nrow(trajA1))
f03 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f03 |> println
f13 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f13 |> println

β_ = range(-2, 3, length = 2001)
M0 = [lorenz_matrix(σ_[k], ρ_[k], b_[k]) for k in eachindex(β_)]
M1 = [syntheticSINDy((1-β)*f01.matrix + β*f11.matrix, vrbl, cnfg, method = "SINDy").matrix for β in β_]
M2 = [syntheticSINDy((1-β)*f02.matrix + β*f12.matrix, vrbl, cnfg, method = "SINDy").matrix for β in β_]
M3 = [syntheticSINDy((1-β)*f03.matrix + β*f13.matrix, vrbl, cnfg, method = "SINDy").matrix for β in β_]

E1 = rmse.(M0, M1)
E2 = rmse.(M0, M2)
E3 = rmse.(M0, M3)

plot(β_, E1, legend = :none, yscale = :log10, color = :black, xticks = [βm, β0, β1, βM], xlims = [βm, βM], size = [300, 150], ylims = [1e-15, 1e-12], yticks = [1e-15, 1e-12]); png("temp5")

plot(size = [600, 150], yscale = :log10, ylims = [1e-3, 1e+2], yticks = [1e-3, 1e+0, 1e+2], left_margin = -0.2mm)
plot!(β_[1:2:end], E2[1:2:end], color = 3, ms = 1, msw = 0, lw = 2)
plot!(β_[1:2:end], E3[1:2:end], color = 4, ms = 1, msw = 0, lw = 2)
png("temp")

E2 |> maximum
E3 |> minimum

png("lorenz matrix error linear")

using MATLAB

mat"""
plot($β_, $E2)
hold on
plot($β_, $E3)
xlabel('beta')
ylabel('E')
legend({'ground truth (7) vs (E3) AFS, with noise η = 10^{-1}', 'ground truth (7) vs (E6) AFS, with noise η = 10^{0}'}, 'Location', 'best')
"""

cnfg = cook(vrbl, poly = 0:2)
ΘX1 = SINDy(_trajA0, vrbl, cnfg; λ = 1e-3).matrix
ΘX2 = SINDy(_trajA1, vrbl, cnfg; λ = 1e-3).matrix
η_l = range(-3, 2, length = 1001)
norm_μ1 = []
norm_μ2 = []
@showprogress for η = η_l
    Random.seed!(1)
    saveat = 1900:1e-3:2000
    trajA0 = deepcopy(_trajA0); trajA0.x .+= exp10(η)*randn(nrow(trajA0)); trajA0.y .+= exp10(η)*randn(nrow(trajA0)); trajA0.z .+= exp10(η)*randn(nrow(trajA0))
    trajA1 = deepcopy(_trajA1); trajA1.x .+= exp10(η)*randn(nrow(trajA1)); trajA1.y .+= exp10(η)*randn(nrow(trajA1)); trajA1.z .+= exp10(η)*randn(nrow(trajA1))

    f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-3); f0 |> println
    f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-3); f1 |> println
    push!(norm_μ1, rmse(f01.matrix, f0.matrix))
    push!(norm_μ2, rmse(f11.matrix, f1.matrix))
end
plt_η = plot(size = [600, 150], scale = :log10, ylims = [1e-7, 1e+2], yticks = [1e-7, 1e+0, 1e+2], xticks = exp10.(-6:6), left_margin = -0.2mm)#, ylims = [-10, 100], yticks = [0, 100])
plot!(plt_η, exp10.(η_l), norm_μ1, ms = 1, color = 1, lw = 2)
plot!(plt_η, exp10.(η_l), norm_μ2, ms = 1, color = 2, lw = 2)
# plot!(plt_η, exp10.(η_l), norm_μ1, color = 1, xticks = exp10.(-6:6), lw = 2)
# plot!(plt_η, exp10.(η_l), norm_μ2, color = 2, xticks = exp10.(-6:6), lw = 2)
png("temp")


β_ = range(-100, 100, length = 2001)
M0 = [lorenz_matrix(
    σ_[801] + (σ_[1201]-σ_[801])β,
    ρ_[801] + (ρ_[1201]-ρ_[801])β,
    b_[801] + (b_[1201]-b_[801])β
    ) for β in β_]
M1 = [syntheticSINDy((1-β)*f01.matrix + β*f11.matrix, vrbl, cnfg, method = "SINDy").matrix for β in β_]
E1 = rmse.(M0, M1)
plot(size = [300, 150])
plot!(β_, E1, yscale = :log10, ylims = [1e-15, 1e-11], yticks = [1e-15, 1e-11], xticks = [-100, 100], legend = :none, color = :black); png("temp4")
