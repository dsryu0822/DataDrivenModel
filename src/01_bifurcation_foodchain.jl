include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    food chain

''''''''''''''''''''''''''''''''''''''''''''''''''"""
ic1 = [0.85, 0.4, 0.8]
ic2 = [0.85, 0.2, 0.8]
# sol1 = factory_foodchain(DataFrame, 0.9612, ic = ic1, saveat = 9000:1e-1:10000)
# sol2 = factory_foodchain(DataFrame, 0.9612, ic = ic2, saveat = 9000:1e-1:10000)
# plot(
#     plot(sol1.R, sol1.C, sol1.P; color = :blue, lw = 1),
#     plot(sol2.R, sol2.C, sol2.P; color = :blue, lw = 1),
# )
ic_ = [[0.85, 0.3rand() + 0.2, 0.8] for _ in 1:10]

# temp = plot(xlims = [0.960, 0.962])
# for k in findall(0.955 .< K_ .< 0.965)
#     v = bfcn[K_[k]]
#     h = fill(K_[k], length(v))
#     scatter!(temp, h, v)
# end
# temp

pm, pM = .88, 1.00; p0, p1 = 0.93, 0.96;
K_ = range(pm, pM, length = 2001)
bfcn = callbfcn("G:/BF/foodchain/bfcnA.jld2")
@showprogress @threads for k in eachindex(K_)[0.955 .< K_ .< 0.965]
    sol1 = factory_foodchain(DataFrame, K_[k], ic = ic1, saveat = 9000:1e-1:10000)
    sol2 = factory_foodchain(DataFrame, K_[k], ic = ic2, saveat = 9000:1e-1:10000)
    bfcn[K_[k]] = [sol1.P[arglmin(sol1.P)]; sol2.P[arglmin(sol2.P)]]
    if minimum(bfcn[K_[k]]) > 0.55 continue end
    for l in eachindex(ic_)
        ic = ic_[l]
        sol = factory_foodchain(DataFrame, K_[k], ic = ic, saveat = 9000:1e-1:10000)
        P_ = sol.P
        if !isempty(P_) && minimum(P_) > 0.55
            bfcn[K_[k]] = P_[arglmin(P_)]
            break
        end
    end
end
scatter(dict2bifurcation(bfcn)..., xticks = [pm, 0.93, 0.94, 0.95, 0.96, pM], ylims = [0.55, 0.8], yticks = [0.55, 0.8], ms = .3, msw = 0, color = :black, size = [600, 200]); png("temp")
bfcn = cleanse(Float64, bfcn)
JLD2.@save "G:/BF/foodchain/bfcnA.jld2" bfcn

for (P0, P1) = [[93, 96], [94, 96], [95, 96]]
# for (P0, P1) = [[93, 94], [93, 95], [94, 95], [93, 96], [94, 96], [95, 96]]
    # P0 = 95; P1 = 96
    p0, p1 = P0/100, P1/100
    @info "G:/BF/foodchain/bfcnB_$(P0)$(P1).jld2"

    βm = (pm - p0) / (p1 - p0)
    βM = (pM - p0) / (p1 - p0)
    β0, β1 = 0, 1
    β_ = range(βm, βM, length = 2001)

    trajA0 = factory_foodchain(DataFrame, p0, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
    trajA1 = factory_foodchain(DataFrame, p1, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
    vrbl = reverse(half(names(trajA0[:, Not(:t)])))
    
    cnfg = cook(vrbl, poly = 0:4)
    f0 = SINDy(trajA0, vrbl, cnfg; λ = 1e-8); f0 |> println
    f1 = SINDy(trajA1, vrbl, cnfg; λ = 1e-8); f1 |> println
    # trajB0 = ssolve(f0, trajA0[[1], f0.rname], 4000:1e-2:5000)
    # trajB1 = ssolve(f1, trajA1[[1], f1.rname], 4000:1e-2:5000)
    cnfg = cookPI(vrbl, poly = 0:3)
    g0 = SINDyPI(trajA0, vrbl, cnfg; λ = 1e-8); g0 |> println
    g1 = SINDyPI(trajA1, vrbl, cnfg; λ = 1e-8); g1 |> println
    # trajC0 = ssolve(g0, trajA0[[1], g0.rname], 4000:1e-2:5000)
    # trajC1 = ssolve(g1, trajA1[[1], g1.rname], 4000:1e-2:5000)

    # plot(
    #     plot(trajA0.R, trajA0.C, trajA0.P, alpha = .5, color = :black),
    #     plot(trajA1.R, trajA1.C, trajA1.P, alpha = .5, color = :black),
    #     plot(trajB0.R, trajB0.C, trajB0.P, alpha = .5, color = :red),
    #     plot(trajB1.R, trajB1.C, trajB1.P, alpha = .5, color = :red),
    #     plot(trajC0.R, trajC0.C, trajC0.P, alpha = .5, color = :blue),
    #     plot(trajC1.R, trajC1.C, trajC1.P, alpha = .5, color = :blue),
    #     layout = (:, 2), size = [400, 600]
    # )

    f_ = affine(Function, f0, f1)
    bfcn = callbfcn("G:/BF/foodchain/bfcnB_$(P0)$(P1).jld2")
    @showprogress for k in eachindex(β_)[0.955 .< K_ .< 0.965]
        sol1 = solve(ODEProblem(f_, ic1, (0.0, 10000.0), (β_[k],)), RK4(), dt = 1e-1, adaptive=false, maxiters = Inf)
        sol2 = solve(ODEProblem(f_, ic2, (0.0, 10000.0), (β_[k],)), RK4(), dt = 1e-1, adaptive=false, maxiters = Inf)
        if sol1.retcode == ReturnCode.Success && sol2.retcode == ReturnCode.Success
            P1_ = sol1[3, sol1.t .≥ 9000]
            P2_ = sol2[3, sol2.t .≥ 9000]
            bfcn[β_[k]] = [P1_[arglmin(P1_)]; P2_[arglmin(P2_)]]
        end
        if minimum(bfcn[β_[k]]) > 0.10 continue end
        for l in 1:10
            ic = ic_[l]
            sol = solve(ODEProblem(f_, ic, (0.0, 10000.0), (β_[k],)), RK4(), dt = 1e-1, adaptive=false, maxiters = Inf)
            if sol.retcode == ReturnCode.Success
                P_ = sol[3, sol.t .≥ 9000]
                if !isempty(P_) && (minimum(P_) > 0.10)
                    bfcn[β_[k]] = P_[arglmin(P_)]
                    break
                end
            end
        end
    end
    # scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], ms = .5, ma = .5, msw = 0, color = :red, size = [400, 200], xlims = [βm, βM]); png("temp1")
    bfcn = cleanse(Float64, bfcn)
    JLD2.@save "G:/BF/foodchain/bfcnB_$(P0)$(P1).jld2" bfcn

    g_ = affinePI(Function, g0, g1)
    bfcn = callbfcn("G:/BF/foodchain/bfcnC_$(P0)$(P1).jld2")
    @showprogress @threads for k in eachindex(β_)[0.955 .< K_ .< 0.965]
        sol1 = solve(ODEProblem(g_, ic1, (0.0, 10000.0), (β_[k],)), RK4(), dt = 1e-1, adaptive=false, maxiters = Inf)
        sol2 = solve(ODEProblem(g_, ic2, (0.0, 10000.0), (β_[k],)), RK4(), dt = 1e-1, adaptive=false, maxiters = Inf)
        if sol1.retcode == ReturnCode.Success && sol2.retcode == ReturnCode.Success
            P1_ = sol1[3, sol1.t .≥ 9000]
            P2_ = sol2[3, sol2.t .≥ 9000]
            bfcn[β_[k]] = [P1_[arglmin(P1_)]; P2_[arglmin(P2_)]]
        end
        if minimum(bfcn[β_[k]]) > 0.10 continue end
        for l in 1:10
            ic = ic_[l]
            sol = solve(ODEProblem(g_, ic, (0.0, 10000.0), (β_[k],)), RK4(), dt = 1e-1, adaptive=false, maxiters = Inf)
            if sol.retcode == ReturnCode.Success
                P_ = sol[3, sol.t .≥ 9000]
                if !isempty(P_) && (minimum(P_) > 0.10)
                    bfcn[β_[k]] = P_[arglmin(P_)]
                    break
                end
            end
        end
    end
    # # scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .5, msw = 0, color = :blue); png("temp2")
    bfcn = cleanse(Float64, bfcn)
    JLD2.@save "G:/BF/foodchain/bfcnC_$(P0)$(P1).jld2" bfcn
end

"""''''''''''''''''''''''''''''''''''''''''''''''''''

            pitchfork bifurcation

''''''''''''''''''''''''''''''''''''''''''''''''''"""
sol = factory_foodchain(DataFrame, 0.9612, ic = [.85, .4, .8], saveat = 9000:1e-1:9500)
ic1 = collect(sol[argmin(sol.P), [:R, :C, :P]])
sol = factory_foodchain(DataFrame, 0.9612, ic = [.85, .2, .8], saveat = 9000:1e-1:9500)
ic2 = collect(sol[argmin(sol.P), [:R, :C, :P]])

# K_ = 0.96080:1e-5:0.96150
K_ = [0.96, 0.960802, 0.961, 0.961531, 0.962]
plt_ = []
@showprogress for k in eachindex(K_)
    sol = factory_foodchain(DataFrame, K_[k], ic = ic1, saveat = 0:1e-1:3000)[(end-1000):end, :]
    plt = plot(legend = :none, camera = [45, 15], formatter = _ -> "",
        xlims = extrema(sol.R) .+ (-0.04, 0.04), xticks = [extrema(sol.R)...] .+ (-0.04, 0.04),
        ylims = extrema(sol.C) .+ (-0.04, 0.04), yticks = [extrema(sol.C)...] .+ (-0.04, 0.04),
        zlims = extrema(sol.P) .+ (-0.04, 0.04), zticks = [extrema(sol.P)...] .+ (-0.04, 0.04),
    )
    CSV.write("G:/BF/foodchain/trajA_1$(k).csv", sol)
    if k != 5
        plot!(plt, sol.R, sol.C, sol.P; color = :blue, lw = 1)
    end
    sol = factory_foodchain(DataFrame, K_[k], ic = ic2, saveat = 0:1e-1:3000)[(end-1000):end, :]
    CSV.write("G:/BF/foodchain/trajA_2$(k).csv", sol)
    if k != 1
        plot!(plt, sol.R, sol.C, sol.P; color = :red, lw = 1)
    end
    push!(plt_, plt)
    # png(plt, "G:/K_$(rpad(K_[k], 8, '0')).png")
end
plot(plt_..., layout = (1, :), size = (800, 150), margin = - 0mm)
png("temp")
plot(plt_[1], size = [400, 400], margin = - 0mm)

pm, pM = .9606, 0.9616; p0, p1 = .950, .960;
K_ = range(pm, pM, length = 1001)
bfcn1 = callbfcn()
bfcn2 = callbfcn()
@showprogress @threads for k in eachindex(K_)
    sol1 = factory_foodchain(DataFrame, K_[k], ic = ic1, saveat = 0:1e-1:10000)
    sol2 = factory_foodchain(DataFrame, K_[k], ic = ic2, saveat = 0:1e-1:10000)
    z1_ = sol1.P[sol1.t .≥ 9000]
    z2_ = sol2.P[sol2.t .≥ 9000]
    bfcn1[K_[k]] = z1_[arglmin(z1_)]
    bfcn2[K_[k]] = z2_[arglmin(z2_)]
end
# # JLD2.@save "G:/BF/foodchain/bfcnA_1_9606_9616.jld2" bfcn1
# # JLD2.@save "G:/BF/foodchain/bfcnA_2_9606_9616.jld2" bfcn2
bfcn1 = JLD2.load("G:/BF/foodchain/bfcnA_1_9606_9616.jld2")["bfcn1"]
bfcn2 = JLD2.load("G:/BF/foodchain/bfcnA_2_9606_9616.jld2")["bfcn2"]
for (k, v) in bfcn1
    if k < 0.9614812
        _bfcn1[k] = v
    end
end
plot(size = [100, 100], ticks = [], ylims = [0.55, 0.8])
scatter!(dict2bifurcation(bfcn2)..., ms = .5, ma = .5, msw = 0, color = :red);
scatter!(dict2bifurcation(1bfcn1)..., ms = .5, ma = .5, msw = 0, color = :blue); png("temp")
istaskdone(task)


"""''''''''''''''''''''''''''''''''''''''''''''''''''

            different initial conditions

''''''''''''''''''''''''''''''''''''''''''''''''''"""

# ic_ = [collect(sol[k, [:R, :C, :P]]) for k in shuffle(1:nrow(sol))[1:100]]
Random.seed!(2)
ic_ = DataFrame(R = rand(10), C = rand(10), P = rand(10))
mat"""
nIDs = 26;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}


fig = tiledlayout(5, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
set(gcf, 'Position', [100, 100, 750, 900])
"""

K_ = [0.93, 0.94, 0.95, 0.96]
for i = 1:4
    trajA_, trajB_, trajC_ = [], [], []

    trajA = factory_foodchain(DataFrame, K_[i], ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
    push!(trajA_, trajA)
    push!(trajA_, factory_foodchain(DataFrame, K_[i], ic = [0.820915, 0.158239, 0.01], saveat = 900:1e-1:1000))
    vrbl = reverse(half(names(trajA[:, Not(:t)])))

    cnfg = cook(vrbl, poly = 0:4)
    f0 = SINDy(trajA, vrbl, cnfg; λ = 1e-8)
    cnfg = cookPI(vrbl, poly = 0:3)
    g0 = SINDyPI(trajA, vrbl, cnfg; λ = 1e-8)

    @showprogress for (k, ic) = enumerate(eachrow(ic_))
        trajB = ssolve(f0, ic, 1000:1e-1:3000)
        trajC = ssolve(g0, ic, 1000:1e-1:3000)
        # push!(trajA_, trajA)
        push!(trajB_, trajB)
        push!(trajC_, trajC)
    end
    
    MtrajAR_ = [trajA.R for trajA in trajA_]
    MtrajAC_ = [trajA.C for trajA in trajA_]
    MtrajAP_ = [trajA.P for trajA in trajA_]

    MtrajBR_ = [trajB.R for trajB in trajB_]
    MtrajBC_ = [trajB.C for trajB in trajB_]
    MtrajBP_ = [trajB.P for trajB in trajB_]

    MtrajCR_ = [trajC.R for trajC in trajC_]
    MtrajCC_ = [trajC.C for trajC in trajC_]
    MtrajCP_ = [trajC.P for trajC in trajC_]

    @mput MtrajAR_ MtrajAC_ MtrajAP_ MtrajBR_ MtrajBC_ MtrajBP_ MtrajCR_ MtrajCC_ MtrajCP_
    mat"""
    nexttile($i)
    hold on;
    grid on;
    for k = 1:numel(MtrajAR_)
        plot3(MtrajAR_{k}, MtrajAC_{k}, MtrajAP_{k}, 'Color', [0, 0, 0], 'LineWidth', 0.5);
    end
    xlabel('R'); ylabel('C'); zlabel('P'); xlim([0 0.95]); ylim([0 0.95]); zlim([0 1.2]);
    xticks([0 0.95]); yticks([0 0.95]); zticks([0 1.2]);
    xticklabels([]); yticklabels([]); zticklabels([]);
    ax = gca;
    ax.XColor = 'b';
    view(15, 30)
    if $i == 1
        text(-0.1, 1.0, 0.2, charlbl{$i}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
    end

    nexttile($i + 4)
    hold on;
    grid on;
    for k = 1:numel(MtrajBR_)
        plot3(MtrajBR_{k}, MtrajBC_{k}, MtrajBP_{k}, 'Color', [1, 0, 0], 'LineWidth', 0.5);
    end
    xlabel('R'); ylabel('C'); zlabel('P'); xlim([0 0.95]); ylim([0 0.95]); zlim([0 1.2]);
    xticks([0 0.95]); yticks([0 0.95]); zticks([0 1.2]);
    xticklabels([]); yticklabels([]); zticklabels([]);
    view(15, 30)
    if $i == 1
        text(-0.1, 1.0, 0.2, charlbl{$i+1}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
    end

    nexttile($i + 8)
    hold on;
    grid on;
    for k = 1:numel(MtrajAR_)
        plot3(MtrajAR_{k}, MtrajAC_{k}, MtrajAP_{k}, 'Color', [0, 0, 0], 'LineWidth', 0.5);
    end
    for k = 1:numel(MtrajBR_)
        plot3(MtrajBR_{k}, MtrajBC_{k}, MtrajBP_{k}, 'Color', [1, 0, 0], 'LineWidth', 0.5);
    end
    xlabel('R'); ylabel('C'); zlabel('P'); xlim([0 0.95]); ylim([0 0.95]); zlim([0 1.2]);
    xticks([0 0.95]); yticks([0 0.95]); zticks([0 1.2]);
    xticklabels([]); yticklabels([]); zticklabels([]);
    view(15, 30)
    if $i == 1
        text(-0.1, 1.0, 0.2, charlbl{$i+2}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
    end

    nexttile($i + 12)
    hold on;
    grid on;
    for k = 1:numel(MtrajCR_)
        plot3(MtrajCR_{k}, MtrajCC_{k}, MtrajCP_{k}, 'Color', [1, 0, 0], 'LineWidth', 0.5);
    end
    xlabel('R'); ylabel('C'); zlabel('P'); xlim([0 0.95]); ylim([0 0.95]); zlim([0 1.2]);
    xticks([0 0.95]); yticks([0 0.95]); zticks([0 1.2]);
    xticklabels([]); yticklabels([]); zticklabels([]);
    view(15, 30)
    if $i == 1
        text(-0.1, 1.0, 0.2, charlbl{$i+3}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
    end

    nexttile($i + 16)
    hold on;
    grid on;
    for k = 1:numel(MtrajAR_)
        plot3(MtrajAR_{k}, MtrajAC_{k}, MtrajAP_{k}, 'Color', [0, 0, 0], 'LineWidth', 0.5);
    end
    for k = 1:numel(MtrajCR_)
        plot3(MtrajCR_{k}, MtrajCC_{k}, MtrajCP_{k}, 'Color', [1, 0, 0], 'LineWidth', 0.5);
    end
    xlabel('R'); ylabel('C'); zlabel('P'); xlim([0 0.95]); ylim([0 0.95]); zlim([0 1.2]);
    xticks([0 0.95]); yticks([0 0.95]); zticks([0 1.2]);
    xticklabels([]); yticklabels([]); zticklabels([]);
    view(15, 30)
    if $i == 1
        text(-0.1, 1.0, 0.2, charlbl{$i+4}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
    end
    """
end


labels = ["($a)" for a in 'a':'l']
plot([plot(rand(10)) for _ in eachindex(labels)]..., title = reshape(labels, 1, :), titleloc = :left, titlefont = font(10))
# plot([plt_attA_; plt_attB_; plt_attC_]...)
plot(plt_attA, plt_attB, plt_attC; layout = (1, 3), size = (1200, 400)); png("temp")

trajA_ = [factory_foodchain(DataFrame, p, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000) for p in [0.93, 0.94, 0.95, 0.96]]
plot(
    [plot(df.C, df.R, df.P;
    color = :black, formatter = _ -> "", camera = [75, 30],
    xticks = [minimum(df.C)-0.01], yticks = [maximum(df.R)+0.013], zticks = [minimum(df.P)-0.01]
    ) for df in trajA_]...,
    layout = (1, 4), size = (800, 200)
); png("temp")


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    model reduction

''''''''''''''''''''''''''''''''''''''''''''''''''"""


for λ in [-4, -3, -2, -1]
# for λ in [-8]
    bfcnA = JLD2.load("G:/BF/foodchain/bfcnA.jld2")["bfcn"]
    
    plt_ = []
    for (P0, P1) = [[93, 94], [93, 95], [94, 95], [93, 96], [94, 96], [95, 96]]
        p0, p1 = P0/100, P1/100
        @info "G:/BF/foodchain/bfcnB_$(P0)$(P1)_$(λ).jld2"
        bfcn = JLD2.load("G:/BF/foodchain/bfcnB_$(P0)$(P1)_$(λ).jld2")["bfcn"]

        pm = 0.88; pM = 1.00;
        βm = (pm - p0) / (p1 - p0)
        βM = (pM - p0) / (p1 - p0)
        
        plt = scatter(dict2bifurcation(bfcnA)..., ms = .5, ma = .5, msw = 0, color = :black, xlims = [0.88, 1.00], ylims = [0.55, 0.8], yticks = [0.55, 0.8], xticks = [0.88, p0, p1, 1.00])
        scatter!(twiny(plt), dict2bifurcation(bfcn)..., xticks = [], ms = .5, ma = .5, msw = 0, color = :red, ylims = [0.55, 0.8], xlims = [βm, βM], yticks = [0.55, 0.8])
        push!(plt_, plt)
    end
    plot(plt_..., layout = (2, :), size = [1200, 600]);
    png("SINDy_$(λ).png")
end