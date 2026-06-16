include("../core/header.jl")
include("../core/factorio.jl")

sol = factory_lorenz63(DataFrame, [10, 28, 8/3])[200000:end, :]
plot(sol.x, sol.y, sol.z, alpha = .5)

σ_ = range(6, 15, length = 1001)
ρ_ = range(120, 150, length = 1001)
β_ = range(3, 5, length = 1001)
if !isfile("bifurcation_lorenz63.jld2")
    @info "Calculating bifurcation data for Lorenz63..."
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(σ_)
        sol = factory_lorenz63(DataFrame, [σ_[k], ρ_[k], β_[k]], ic = [100, 100, 100])
        z_ = sol.z[sol.t .≥ 900]
        bfcn[k] = z_[arglmax(z_)]
    end
    JLD2.@save "bifurcation_lorenz63.jld2" bfcn
else
    @info "Loading bifurcation data for Lorenz63 from file..."
    JLD2.@load "bifurcation_lorenz63.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .1, msw = 0, color = :black); png("temp")

traj0 = factory_lorenz63(DataFrame, [σ_[1], ρ_[1], β_[1]], saveat = 900:1e-3:1000)
traj1 = factory_lorenz63(DataFrame, [σ_[101], ρ_[101], β_[101]], saveat = 900:1e-3:1000)

vrbl = reverse(half(names(traj0[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:2)
f0 = SINDy(traj0, vrbl, cnfg; λ = 1e-3); f0 |> print
f1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); f1 |> print

rcvd0 = ssolve(f0, rand(3), 900:1e-3:1000)
rcvd1 = ssolve(f1, rand(3), 900:1e-3:1000)
plot(
    plot(traj0.x, traj0.y, traj0.z, color = :black),
    plot(traj1.x, traj1.y, traj1.z, color = :black),
    plot(rcvd0[:, 2], rcvd0[:, 3], rcvd0[:, 4], color = :blue),
    plot(rcvd1[:, 2], rcvd1[:, 3], rcvd1[:, 4], color = :blue),
)

if !isfile("bifurcation_lorenz63_recover.jld2")
    @info "Calculating bifurcation data for Lorenz63 from recovered models..."
    α_ = 0:1e-2:10
    f_ = [syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg, method = "SINDy") for α in α_]
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(α_)
        sol = solve(ODEProblem(define(Function, f_[k]; fname = "f_$k"), [100, 100, 100], (0, 1000)), RK4(), maxiters = 1e+7, dt = 1e-3, adaptive=false, saveat = 900:1e-3:1000)
        z_ = sol[3, :]
        bfcn[α_[k]] = z_[arglmax(z_)]
    end
    JLD2.@save "bifurcation_lorenz63_recover.jld2" bfcn
else
    @info "Loading bifurcation data for Lorenz63 from recovered models from file..."
    JLD2.@load "bifurcation_lorenz63_recover.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .1, msw = 0, color = :blue); png("temp")

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                        food chain

''''''''''''''''''''''''''''''''''''''''''''''''''"""
sol = factory_foodchain(DataFrame, 0.98, saveat = 0:1e-2:1000)
plot(sol.R, sol.C, sol.P)

K_ = 0.90:1e-4:1.0
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
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, color = :black)
png("temp")

largs = (; xlabel = L"R", ylabel = L"C", zlabel = L"P", xticks = [0.3, 0.8], yticks = [0.2, 0.45], zticks = [0.65, 1.0])

traj0 = factory_foodchain(DataFrame, 0.95, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
traj1 = factory_foodchain(DataFrame, 0.96, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000)
plt_pp_1 = plot(traj0.R, traj0.C, traj0.P; largs..., color = :black)
plt_pp_2 = plot(traj1.R, traj1.C, traj1.P; largs..., color = :black)

vrbl = reverse(half(names(traj0)[2:end]))
cnfg = cook(vrbl; poly = 0:4)

@time f0 = SINDy(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print
@time f1 = SINDy(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print

tspan = (0, 5000)
prob0 = ODEProblem(define(Function, f0), collect(traj0[1, vrbl[2]]), tspan, RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false; saveat = 4000:1e-2:5000)
sol0 = solve(prob0);
prob1 = ODEProblem(define(Function, f1), collect(traj1[1, vrbl[2]]), tspan, RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false; saveat = 4000:1e-2:5000)
sol1 = solve(prob1);
plt_pp_3 = plot(sol0[1,:], sol0[2, :], sol0[3, :]; largs..., color = :red)
plt_pp_4 = plot(sol1[1,:], sol1[2, :], sol1[3, :]; largs..., color = :red)
plot(plt_pp_1, plt_pp_2, plt_pp_3, plt_pp_4)

if !isfile("bifurcation_foodchain_failed.jld2")
    @info "Calculating bifurcation data for foodchain from recovered models..."
    α_ = -6:1e-2:5
    f_ = [define(Function, syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg, method = "SINDy"), fname = "f_$k") for (k, α) in enumerate(α_)]
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(α_)        
        for _ in 1:1000
            try
            sol = solve(ODEProblem(f_[k], collect(traj1[1, vrbl[2]]), (0, 10000), RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false; saveat = 0:1e-1:10000))
            P_ = sol[3, sol.t .≥ 9000]
            if !isempty(P_) && minimum(P_) > 0.55
                bfcn[α_[k]] = P_[arglmin(P_)]
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

vrbl = reverse(half(names(traj0)[2:end]))
cnfg = cookPI(vrbl; poly = 0:3)
@time f0 = SINDyPI(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print # define(f0) |> print
@time f1 = SINDyPI(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print # define(f1) |> print

prob0 = DAEProblem(define(Function, f0), collect(traj0[1, vrbl[1]]), collect(traj0[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1])))
sol0 = solve(prob0, Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6);
prob1 = DAEProblem(define(Function, f1), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1])))
sol1 = solve(prob1, Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6)
plt_pp_5 = plot(eachrow(stack(sol0.u))...; color = :blue)
plt_pp_6 = plot(eachrow(stack(sol1.u))...; color = :blue)

if !isfile("bifurcation_foodchain_recover.jld2")
    @info "Calculating bifurcation data for foodchain from recovered models..."
    α_ = -6:1e-2:5
    f_ = [define(Function, syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg, method = "SINDyPI"), fname = "f_$k") for (k, α) in enumerate(α_)]
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(α_)
        for _ in 1:1000
            sol = solve(DAEProblem(f_[k], collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), (0, 10000), differential_vars = ones(Bool, length(vrbl[1]))), Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6; saveat = 0:1e-1:10000)
            P_ = sol[3, sol.t .≥ 9000]
            if !isempty(P_) && minimum(P_) > 0.55
                bfcn[α_[k]] = P_[arglmin(P_)]
                break
            end
        end
    end
    JLD2.@save "bifurcation_foodchain_recover.jld2" bfcn
else
    @info "Loading bifurcation data for foodchain from recovered models from file..."
    JLD2.@load "bifurcation_foodchain_recover.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = 1, ma = .5, msw = 0, color = :blue); png("temp")

_bfcn = deepcopy(bfcn)

bfcn = deepcopy(_bfcn)