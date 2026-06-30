include("../core/header.jl")
include("../core/factorio.jl")
condition(u,t,integrator) = 0.1 - u[3]
affect!(integrator) = terminate!(integrator)
early_stop = ContinuousCallback(condition,affect!)
solveargs = (; reltol = 1e-24, initializealg = DiffEqBase.BrownFullBasicInit(), callback = early_stop)
function foodchain(du, u, p, t)
    R,C,P = u
    K = p[1]
     xc,    yc,   xp,    yp,      R0,  C0 = (
    0.4, 2.009, 0.08, 2.876, 0.16129, 0.5)
    du[1] = R*(1 - frac(R,K)) - xc*yc*frac(C*R, R + R0)
    du[2] = xc*C*(frac(yc*R, R + R0) - 1) - xp*yp*frac(P*C, C + C0)
    du[3] = xp*P*(frac(yp*C, C + C0) - 1)
    return du
end

traj0[1, 2:4] |> collect
ic_chaos = [.85, .2+.1rand(), .85]
ic_chaos = [0.85, 0.21202661351046087, 0.85]
traj0 = factory_foodchain(DataFrame, 0.99976 - 4e-5, ic = ic_chaos, saveat = 0:1e-1:5000)[1:10:end,:]
traj1 = factory_foodchain(DataFrame, 0.99976 + 4e-5, ic = ic_chaos, saveat = 0:1e-1:5000)[1:10:end,:]
f0 = SINDyPI(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print # define(f0) |> print
f1 = SINDyPI(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print # define(f1) |> print
prob0 = DAEProblem(define(Function, f0), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1])))
sol0 = solve(prob0, Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6);
prob1 = DAEProblem(define(Function, f1), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1])))
sol1 = solve(prob1, Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6)
plot(
    plot(traj0.P, color = :black, ylims = [0, 1.1]),
    plot(traj1.P, color = :black, ylims = [0, 1.1]),
    plot(sol0.t, sol0[3, :], color = :blue, ylims = [0, 1.1]),
    plot(sol1.t, sol1[3, :], color = :blue, ylims = [0, 1.1]),
    xticks = [0, 5000], yticks = [0, .55, 1.1]
)
icc = shuffle([[r...] for r in eachrow(traj0[:, [:R, :C, :P]])])

ε_ = exp10.(range(-5, -2, 20))
tau1__ = [zeros(300) for k in eachindex(ε_)]
for k in eachindex(ε_)
    @showprogress @threads for l in eachindex(tau1__[k])
        prob = ODEProblem(foodchain, icc[l], (0, 100000), [0.99976+ε_[k]])
        sol = solve(prob; solveargs...)
        tau = findfirst(sol[3, :] .< 0.55)
        tau1__[k][l] = !isnothing(tau) ? sol.t[tau] : last(sol.t)
    end
end
scatter(ε_, mean.(tau1__), scale = :log10, size = [400, 400], color = :white, legend = :none, xticks = exp10.(-10:-1))

traj0 = factory_foodchain(DataFrame, 0.95, ic = [0.820915, 0.158239, 0.953786], saveat = 0:1e-1:5000)
traj1 = factory_foodchain(DataFrame, 0.96, ic = [0.820915, 0.158239, 0.953786], saveat = 0:1e-1:5000)
vrbl = reverse(half(names(traj0)[2:end]))
cnfg = cookPI(vrbl; poly = 0:3)
f0 = SINDyPI(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print
f1 = SINDyPI(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print

α_c = 4.777865647
icc = shuffle(solve(DAEProblem(define(Function, syntheticSINDy((1-α_c)*f0.matrix + α_c*f1.matrix, vrbl, cnfg[1], method = "SINDyPI")), ones(3), [0.820915, 0.158239, 0.953786], (0, 1e+5), differential_vars = ones(Bool, length(vrbl[1]))), Sundials.IDA(); solveargs...).u)

α_ = α_c .+ 40ε_
f_ = []
for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg[1], method = "SINDyPI")
    push!(f_, define(Function, g; fname = "f_$k"))
end
tau2__ = [zeros(500) for k in eachindex(α_)]
for k in eachindex(α_)
    @showprogress @threads for l in eachindex(tau2__[k])
        probg = DAEProblem(f_[k], ones(3), icc[l], (0, 1e+5), differential_vars = ones(Bool, length(vrbl[1])))
        solg = solve(probg, Sundials.IDA(); solveargs...)
        t_ = solg.t
        u_ = solg.u
        for _ in 1:100
            if last(u_)[3] > .55
                probg = DAEProblem(f_[k], ones(3), icc[l], (last(t_), 1e+5), differential_vars = ones(Bool, length(vrbl[1])))
                solg = solve(probg, Sundials.IDA(); solveargs...)
                t_ = solg.t; u_ = solg.u
            end
        end
        tau = findfirst(stack(u_)[3, :] .< 0.55)
        tau2__[k][l] = !isnothing(tau) ? t_[tau] : last(t_)
    end
end
# @JLD2.save "transient.jld2" ε_ tau1__ tau2__
# @JLD2.load "transient.jld2"
plt_trs = plot(size = [400, 400], xticks = exp10.(-10:-1), ylims = [1e+2, 1e+5]);
scatter!(plt_trs, ε_, mean.(tau1__), scale = :log10, color = :white, msc = :black, shape = :circ,  label = "actual");
scatter!(plt_trs, ε_, mean.(tau2__), scale = :log10, color = :white, msc = :blue,  shape = :rect, label = "recovered");
plt_trs
png("temp")

maximum.(tau2__)