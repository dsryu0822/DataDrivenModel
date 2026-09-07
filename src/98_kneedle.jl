using Kneedle

x = 1:100
y = cumsum(1 ./ rand(100))

kr = kneedle(x, y)
knees(kr)
plot(x, y)
vline!(knees(kneedle(x, y, "concave_inc", 1)))

eachcol(f1.matrix)

f0 = SINDy(trajA1, vrbl, cnfg)
Matrix(f0.matrix)
f0 = SINDy(trajA1, vrbl, cnfg, kneedle = true; λ = 1e-1)
Matrix(f0.matrix)
Matrix(f0.matrixF)

pm, pM = .88, 1.00;
    P0 = 95; P1 = 96
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

    λ = -8
    f0 = SINDy(trajA0, vrbl, cnfg, kneedle = true; λ = exp10(λ)); f0 |> println
    f1 = SINDy(trajA1, vrbl, cnfg, kneedle = true; λ = exp10(λ)); f1 |> println

    trajB0 = ssolve(f0, trajA0[[1], f0.rname], 4000:1e-2:5000)
    trajB1 = ssolve(f1, trajA1[[1], f1.rname], 4000:1e-2:5000)

    plot(
        plot(trajA0.R, trajA0.C, trajA0.P, alpha = .5, color = :black),
        plot(trajA1.R, trajA1.C, trajA1.P, alpha = .5, color = :black),
        plot(trajB0.R, trajB0.C, trajB0.P, alpha = .5, color = :red),
        plot(trajB1.R, trajB1.C, trajB1.P, alpha = .5, color = :red),
        layout = (:, 2), size = [400, 600]
    )

    
    cnfg = cook(vrbl, poly = 0:4)
    f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
    bfcn = callbfcn()
    @showprogress @threads for k in eachindex(β_)
        for l in 1:10
            sol = ssolve(f_[k], trajA0[[1+100(l-1)], f_[k].rname], 0:1e-1:10000)
            P_ = sol.P[sol.t .≥ 9000]
            if !isempty(P_) && (minimum(P_) > 0.10)
                bfcn[β_[k]] = P_[arglmin(P_)]
                break
            end
        end
    end
    # scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], ms = .5, ma = .5, msw = 0, color = :red, size = [400, 200], xlims = [βm, βM]); png("temp1")
    _bfcn = deepcopy(bfcn)
    bfcn = deepcopy(_bfcn)
    JLD2.@save "G:/BF/foodchain/bfcnB_$(P0)$(P1)_$(λ).jld2" bfcn