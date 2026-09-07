include.("../core/" .* readdir("core")[[1,2,3,4,6]])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Lorenz-63

''''''''''''''''''''''''''''''''''''''''''''''''''"""
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

define(f0) |> print
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
    β_ = 0:1e-2:10
    f_ = [syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy") for β in β_]
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(β_)
        sol = solve(ODEProblem(define(Function, f_[k]; fname = "f_$k"), [100, 100, 100], (0, 1000)), RK4(), maxiters = 1e+7, dt = 1e-3, adaptive=false, saveat = 900:1e-3:1000)
        z_ = sol[3, :]
        bfcn[β_[k]] = z_[arglmax(z_)]
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
# sol = factory_foodchain(DataFrame, 1, ic = [.85, .14873587168078906, .8], saveat = 0:1e-1:10000)
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

begin
    R0 = .16129; Kc = .99976; K0 = .95; K1 = .96
    fβ(K0, K1) = (frac(1,K0) - frac(1,Kc)) / (frac(1,K0) - frac(1,K1))
    K0K1β = []
    for K0 in 0.9:1e-2:1.0
        for K1 in 0.9:1e-2:1.0
            if K0 ≥ K1 continue end
            βmin = (0.9 - K0) / (K1 - K0)
            βmax = (1.0 - K0) / (K1 - K0)
            push!(K0K1β, (K0, K1, fβ(K0, K1), βmin, βmax))
        end
    end
    apc = DataFrame(stack(K0K1β)', [:K0, :K1, :β, :βmin, :βmax])
    apc.βrel = (apc.β - apc.βmin) ./ (apc.βmax - apc.βmin)
    # CSV.write("apc.csv", apc, bom = true)
end

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                        Thomas

''''''''''''''''''''''''''''''''''''''''''''''''''"""
ic0 = [0.0, 0.0, .1]

sol = factory_thomas(DataFrame, 0.1)
plot(sol.x, sol.y, sol.z)

b_ = .14:2e-5:.16
if !isfile("bifurcation_thomas.jld2")
    @info "Calculating bifurcation data for Thomas..."
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(b_)
        sol = factory_thomas(DataFrame, b_[k], ic = ic0)
        x_ = sol.x[sol.t .≥ 1000]
        bfcn[b_[k]] = x_[arglmax(x_)]
    end
    JLD2.@save "bifurcation_thomas.jld2" bfcn
else
    @info "Loading bifurcation data for Thomas from file..."
    JLD2.@load "bifurcation_thomas.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, color = :black); png("temp")

traj0 = factory_thomas(DataFrame, .140, saveat = 1000:1e-2:2000)
traj1 = factory_thomas(DataFrame, .141, saveat = 1000:1e-2:2000)
plt_pp_1 = plot(traj0.x, traj0.y, traj0.z, color = :black)
plt_pp_2 = plot(traj1.x, traj1.y, traj1.z, color = :black)

vrbl = reverse(half(names(traj0)[2:end]))
cnfg = cook(vrbl; poly = 0:1, trig = [1], format = sin)

@time f0 = SINDy(traj0, vrbl, cnfg; λ = 1e-3); f0 |> print
@time f1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); f1 |> print

prob0 = ODEProblem(define(Function, f0), collect(traj0[1, vrbl[2]]), (0, 2000))
sol0 = solve(prob0, RK4(), dt = 1e-2, adaptive=false);
plt_pp_3 = plot(eachrow(stack(sol0.u))...; color = :blue)
prob1 = ODEProblem(define(Function, f1), collect(traj1[1, vrbl[2]]), (0, 2000))
sol1 = solve(prob1, RK4(), dt = 1e-2, adaptive=false);
plt_pp_4 = plot(eachrow(stack(sol1.u))...; color = :blue)
plot(plt_pp_1, plt_pp_2, plt_pp_3, plt_pp_4, legend = :none, size = [600, 600])

if !isfile("bifurcation_thomas_recover.jld2")
    @info "Calculating bifurcation data for Thomas from recovered models..."
    β_ = 0:1e-2:20
    f_ = [define(Function, syntheticSINDy((1-β)*f0.matrix + β*f1.matrix, vrbl, cnfg, method = "SINDy"), fname = "f_$k") for (k, β) in enumerate(β_)]
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(β_)
        sol = solve(ODEProblem(f_[k], ic0, (0, 2000)), RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false)
        x_ = sol[1, sol.t .> 1000]
        bfcn[β_[k]] = x_[arglmax(x_)]
    end
    JLD2.@save "bifurcation_thomas_recover.jld2" bfcn
else
    @info "Loading bifurcation data for Thomas from recovered models from file..."
    JLD2.@load "bifurcation_thomas_recover.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, color = :blue); png("temp")

# ic0 = [0.0, 0.0, .1]
# traj0 = factory_thomas(DataFrame, .135, ic = ic0, saveat = 1000:1e-2:2000)
cnfg = cook(vrbl; poly = 0:5)
traj0 = factory_thomas(DataFrame, .140, ic = ic0, saveat = 1000:1e-2:2000)
traj1 = factory_thomas(DataFrame, .141, ic = ic0, saveat = 1000:1e-2:2000)
@time g0 = SINDy(traj0, vrbl, cnfg; λ = 1e-3); g0 |> print
@time g1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); g1 |> print
plt_pp_1 = plot(traj0.x, traj0.y, traj0.z, color = :black)
plt_pp_2 = plot(traj1.x, traj1.y, traj1.z, color = :black)

if !isfile("bifurcation_thomas_failed.jld2")
    @info "Calculating bifurcation data for Thomas from recovered models with failed SINDy..."
    β_ = 0:1e-2:20
    g_ = [define(Function, syntheticSINDy((1-β)*g0.matrix + β*g1.matrix, vrbl, cnfg, method = "SINDy"), fname = "g_$k") for (k, β) in enumerate(β_)]
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(β_)
        sol = solve(ODEProblem(g_[k], ic0, (0, 2000)), RK4(), maxiters = 1e+7, dt = 1e-2, adaptive=false)
        x_ = sol[1, sol.t .> 1000]
        bfcn[β_[k]] = x_[arglmax(x_)]
    end
    JLD2.@save "bifurcation_thomas_failed.jld2" bfcn
else
    @info "Loading bifurcation data for Thomas from recovered models with failed SINDy from file..."
    JLD2.@load "bifurcation_thomas_failed.jld2" bfcn
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .5, msw = 0, color = :red); png("temp")

 - f0.matrix + f1.matrix
import Base: show

function show(io::IO, f::Float64)
    @printf(io, "%.16f", f)
end

f0.recipe
maximum(keys(bfcn))
