include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    food chain

''''''''''''''''''''''''''''''''''''''''''''''''''"""
sol = factory_foodchain(DataFrame, .99, ic = [.85, .2, .8], saveat = 0:1e-1:10000)
plot(sol.R, sol.C, sol.P)

pm, pM = .88, 1.00; p0, p1 = .950, .960;
K_ = range(pm, pM, length = 1001)
bfcn = callbfcn("G:/BF/foodchain/bfcnA.jld2")
@showprogress @threads for k in eachindex(K_)
    sol = factory_foodchain(DataFrame, K_[k], ic = [.85, .6, .8], saveat = 0:1e-1:10000)
    z_ = sol.P[sol.t .≥ 9000]
    bfcn[K_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [pm, p0, p1, pM], ms = .5, ma = .5, msw = 0, color = :black); png("temp")
# JLD2.@save "G:/BF/foodchain/bfcnA.jld2" bfcn

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

plot(
    plot(trajA0.R, trajA0.C, trajA0.P, alpha = .5, color = :black),
    plot(trajA1.R, trajA1.C, trajA1.P, alpha = .5, color = :black),
    plot(trajB0.R, trajB0.C, trajB0.P, alpha = .5, color = :red),
    plot(trajB1.R, trajB1.C, trajB1.P, alpha = .5, color = :red),
    plot(trajC0.R, trajC0.C, trajC0.P, alpha = .5, color = :blue),
    plot(trajC1.R, trajC1.C, trajC1.P, alpha = .5, color = :blue)
)
CSV.write("G:/BF/foodchain/trajA0.csv", trajA0)
CSV.write("G:/BF/foodchain/trajA1.csv", trajA1)
CSV.write("G:/BF/foodchain/trajB0.csv", trajB0)
CSV.write("G:/BF/foodchain/trajB1.csv", trajB1)
CSV.write("G:/BF/foodchain/trajC0.csv", trajC0)
CSV.write("G:/BF/foodchain/trajC1.csv", trajC1)

βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
β0, β1 = 0, 1
β_ = range(βm, βM, length = 1001)

f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
bfcn = callbfcn("G:/BF/foodchain/bfcnB.jld2")
@showprogress for k in eachindex(β_)
    sol = ssolve(f_[k], trajA0[[1], f_[k].rname], 2000:1e-4:3000)
    z_ = sol.z[sol.t .≥ 2000]
    bfcn[β_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], ms = .5, ma = .5, msw = 0, color = :blue); png("temp")
# JLD2.@save "G:/BF/foodchain/bfcnB.jld2" bfcn

g_ = [syntheticSINDy((1-β)*g0.matrix + β*g1.matrix, vrbl, cnfg, method = "SINDyPI") for β in β_]
bfcn = callbfcn("G:/BF/foodchain/bfcnC.jld2")
@showprogress @threads for k in eachindex(β_)
    sol = ssolve(g_[k], trajA0[[1], g_[k].rname], 900:1e-3:1000)
    z_ = sol.z[sol.t .≥ 900]
    bfcn[β_[k]] = z_[arglmax(z_)]
end
scatter(dict2bifurcation(bfcn)..., xticks = [βm, β0, β1, βM], xlims = [βm, βM], ms = .5, msw = 0, color = :red); png("temp")
# JLD2.@save "G:/BF/foodchain/bfcnC.jld2" bfcn




# ▼ old

sol = factory_foodchain(DataFrame, 1, ic = [.85, .6, .8], saveat = 0:1e-1:10000)
plot(sol.R, sol.C, sol.P)
plot(sol.P)

K_ = 0.88:1e-4:0.90
if !isfile("bifurcation_foodchain.jld2")
    @info "Calculating bifurcation data for food chain..."
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(K_)
        for _ in 1:500
            sol = factory_foodchain(DataFrame, K_[k], saveat = 0:1e-1:10000)
            P_ = sol.P[sol.t .≥ 9000]
            if !isempty(P_) && minimum(P_) > 0.55
                bfcn[K_[k]] = P_[arglmin(P_)]
                break
            end
        end
    end
    JLD2.@save "bifurcation_foodchain.jld2" bfcn
else
    @info "Loading bifurcation data for food chain from file..."
    JLD2.@load "bifurcation_foodchain.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, color = :black, xticks = [.88, p0, p1, 1.0]); png("temp")

largs = (; xlabel = L"R", ylabel = L"C", zlabel = L"P", xticks = [0.3, 0.8], yticks = [0.2, 0.45], zticks = [0.65, 1.0])


p0 = .950; p1 = .960
βmin = (0.88 - p0) / (p1 - p0)
βmax = (1.00 - p0) / (p1 - p0)

traj0 = factory_foodchain(DataFrame, p0, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
traj1 = factory_foodchain(DataFrame, p1, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
plt_pp_1 = plot(traj0.R, traj0.C, traj0.P; largs..., color = :black)
plt_pp_2 = plot(traj1.R, traj1.C, traj1.P; largs..., color = :black)

vrbl = reverse(half(names(traj0)[2:end]))
cnfg = cook(vrbl; poly = 0:4)

@time f0 = SINDy(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print
@time f1 = SINDy(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print

tspan = (0, 5000)
prob0 = ODEProblem(define(Function, f0), collect(traj0[1, vrbl[2]]), tspan, RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false; saveat = 4000:1e-2:5000)
sol0 = solve(prob0, RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false; saveat = 4000:1e-2:5000)
prob1 = ODEProblem(define(Function, f1), collect(traj1[1, vrbl[2]]), tspan, RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false; saveat = 4000:1e-2:5000)
sol1 = solve(prob1, RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false; saveat = 4000:1e-2:5000)
plt_pp_3 = plot(sol0[1,:], sol0[2, :], sol0[3, :]; largs..., color = :red)
plt_pp_4 = plot(sol1[1,:], sol1[2, :], sol1[3, :]; largs..., color = :red)

vrbl = reverse(half(names(traj0)[2:end]))
cnfg = cookPI(vrbl; poly = 0:3); #cnfg[2][[68], :] .= true
@time f0 = SINDyPI(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print # define(f0) |> print
@time f1 = SINDyPI(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print # define(f1, sigdigits = 5) |> print

prob0 = DAEProblem(define(Function, f0), collect(traj0[1, vrbl[1]]), collect(traj0[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1])))
sol0 = solve(prob0, Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6);
prob1 = DAEProblem(define(Function, f1), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1])))
sol1 = solve(prob1, Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6)
plt_pp_5 = plot(eachrow(stack(sol0.u))...; color = :blue)
plt_pp_6 = plot(eachrow(stack(sol1.u))...; color = :blue)
plot(plt_pp_1, plt_pp_2, plt_pp_3, plt_pp_4, plt_pp_5, plt_pp_6; layout = (2, 3), size = (1200, 800))

if !isfile("bifurcation_foodchain_failed.jld2")
    @info "Calculating bifurcation data for foodchain from recovered models..."
    β_ = 2:1e-3:3.78
    f_ = [define(Function, syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy"), fname = "f_$k") for (k, β) in enumerate(β_)]
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(β_)
        for _ in 1:100
            try
            sol = solve(ODEProblem(f_[k], [.85, rand(), .8], (0, 10000), RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false; saveat = 0:1e-1:10000))
            P_ = sol[3, sol.t .≥ 9000]
            if !isempty(P_) && minimum(P_) > 0.55
                bfcn[β_[k]] = P_[arglmin(P_)]
                break
            end
            catch
                continue
            end
        end
    end
    JLD2.@save "bifurcation_foodchain_failed.jld2" bfcn
else
    @info "Loading bifurcation data for foodchain from recovered models from file..."
    JLD2.@load "bifurcation_foodchain_failed.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, color = :red); png("temp")

if !isfile("bifurcation_foodchain_recover.jld2")
    @info "Calculating bifurcation data for foodchain from recovered models..."
    β_ = range(βmin, βmax, length = 1201)
    f_ = [define(Function, syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDyPI"), fname = "f_$k") for (k, β) in enumerate(β_)]
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(β_)
        for _ in 1:100
            sol = solve(DAEProblem(f_[k], zeros(3), [.85, rand(), .8], (0, 10000), differential_vars = ones(Bool, length(vrbl[1]))), Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6; saveat = 0:1e-1:10000)
            P_ = sol[3, sol.t .≥ 9000]
            if !isempty(P_) && minimum(P_) > 0.55
                bfcn[β_[k]] = P_[arglmin(P_)]
                break
            end
        end
    end
    JLD2.@save "bifurcation_foodchain_recover.jld2" bfcn
else
    @info "Loading bifurcation data for foodchain from recovered models from file..."
    JLD2.@load "bifurcation_foodchain_recover.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, color = :blue, xticks = [βmin, 0, 1, βmax]); png("temp")



"""''''''''''''''''''''''''''''''''''''''''''''''''''

            different initial conditions

''''''''''''''''''''''''''''''''''''''''''''''''''"""

# ic_ = [collect(sol[k, [:R, :C, :P]]) for k in shuffle(1:nrow(sol))[1:100]]
ic_ = [rand(3) for _ in 1:100]
plt_attA = plot(legend = :none, lims = [0, 1.2])
plt_attB = plot(legend = :none, lims = [0, 1.2])
plt_attC = plot(legend = :none, lims = [0, 1.2])
@showprogress for ic in ic_
    trajA0 = factory_foodchain(DataFrame, p1, ic = ic, saveat = 900:1e-1:1000)
    trajB0 = ssolve(f1, trajA0[[1], f1.rname], 900:1e-1:1000)
    trajC0 = ssolve(g1, trajA0[[1], g1.rname], 900:1e-1:1000)
    plot!(plt_attA, trajA0.R, trajA0.C, trajA0.P; color = :black)
    plot!(plt_attB, trajB0.R, trajB0.C, trajB0.P; color = :red)
    plot!(plt_attC, trajC0.R, trajC0.C, trajC0.P; color = :blue)
end
plot(plt_attA, plt_attB, plt_attC; layout = (1, 3), size = (1200, 400)); png("temp")