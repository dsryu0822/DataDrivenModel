include("../core/header.jl")
include("tippingutils.jl")

import DifferentialEquations as DE
import Sundials
import DiffEqBase


hrzn = []
vrcl = []
for a0 = 0.15:0.001:2.0
    data = factory_hastingpowell(DataFrame, a0, tspan = 0:1e-2:1000)[50001:end, :]
    bits = arglmax(data.v1)
    push!(hrzn, [a0 for _ in bits])
    push!(vrcl, data.v1[bits])
end
scatter([hrzn...;], [vrcl...;], ms = 1, msw = 0, color = :black, alpha = 0.5, legend = :none, xlims = [1.5, 2.0], ylims = [25, 40])
png("G:/temp.png")

traj0 = factory_hastingpowell(DataFrame, 1.7, tspan = 0:1e-2:1000)[50000:end, :]
traj1 = factory_hastingpowell(DataFrame, 1.8, tspan = 0:1e-2:1000)[50000:end, :]
traj2 = factory_hastingpowell(DataFrame, 2.0, tspan = 0:1e-2:1000)[50000:end, :]
plot(
    plot(traj0.v1, traj0.v2, traj0.v3, color = :black),
    plot(traj1.v1, traj1.v2, traj1.v3, color = :black),
    plot(traj2.v1, traj2.v2, traj2.v3, color = :black),
    legend = :none, layout = (1, 3), size = [1200, 400]
)

vrbl = reverse(half(names(traj0)[Not(1)]))
cnfg = cookPI(vrbl; poly = 0:3)

f0 = SINDyPI(traj0, vrbl, cnfg; λ = 1e-7); f0 |> print
f1 = SINDyPI(traj1, vrbl, cnfg; λ = 1e-7); f1 |> print
f2 = SINDyPI(traj2, vrbl, cnfg; λ = 1e-7); f2 |> print

tspan = (0, 1000)
deargs = (; reltol = 1e-6, initializealg = DiffEqBase.BrownFullBasicInit(), saveat = 0:1e-2:1000)
prob0 = DE.DAEProblem(define(Function, f0), collect(traj0[1, 5:end]), collect(traj0[1, 2:4]), tspan, differential_vars = ones(Bool, length(vrbl[1])))
sol0 = DE.solve(prob0, Sundials.IDA(); deargs...)

prob1 = DE.DAEProblem(define(Function, f1), collect(traj1[1, 5:end]), collect(traj1[1, 2:4]), tspan, differential_vars = ones(Bool, length(vrbl[1])))
sol1 = DE.solve(prob1, Sundials.IDA(); deargs...)

prob2 = DE.DAEProblem(define(Function, f2), collect(traj2[1, 5:end]), collect(traj2[1, 2:4]), tspan , differential_vars = ones(Bool, length(vrbl[1])))
sol2 = DE.solve(prob2, Sundials.IDA(); deargs...)

plot(
    plot(eachrow(stack(sol0.u))..., ticks = [], color = :blue),
    plot(eachrow(stack(sol1.u))..., ticks = [], color = :blue),
    plot(eachrow(stack(sol2.u))..., ticks = [], color = :blue),
    legend = :none, layout = (1, 3), size = [1200, 400]
)

α_ = -2:1e-2:3
hrzn = []
vrtl = []
@showprogress for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg[1], method = "SINDyPI")
    probg = DE.DAEProblem(define(Function, g), collect(traj0[1, 5:end]), collect(traj0[1, 2:4]), (0, 1000), differential_vars = [true, true, true])
    solg = DE.solve(probg, Sundials.IDA(); deargs...)
    traj = stack(solg.u)'[solg.t .≥ 500, :]
    points = arglmax(traj[:, 1])
    push!(vrtl, traj[points, 1])
    push!(hrzn, repeat([α], length(points)))
end
scatter(hrzn, vrtl, ms = 1, msw = 0, color = :blue, alpha = 0.5, legend = :none, xlims = [-2, 3], ylims = [25, 40])
png("temp")



define(f0) |> print
define(f1) |> print

function affined(out, du, u, p, t)
v1, v2, v3 = u; dv1, dv2, dv3 = du
ε = 1e-5
out[1] = (1-ε*t)*(1.7v1 + 0.12v1*v1 - 0.1v1*v2 - 0.005v1*v1*v1 - 0.1v1*dv1 - dv1) + ε*t*(1.8v1 + 0.13v1*v1 - 0.1v1*v2 - 0.005v1*v1*v1 - 0.1v1*dv1 - dv1)
out[2] = (1-ε*t)*(-1.0v2 + 0.1v1*v2 - 0.1v2*v2 - 0.15v2*v3 + 0.01v1*v2*v2 - 0.015v1*v2*v3 - 0.1v1*dv2 - 0.1v2*dv2 - 0.01v1*v2*dv2 - dv2) + ε*t*(-1.0v2 + 0.1v1*v2 - 0.1v2*v2 - 0.15v2*v3 + 0.01v1*v2*v2 - 0.015v1*v2*v3 - 0.1v1*dv2 - 0.1v2*dv2 - 0.01v1*v2*dv2 - dv2)
out[3] = (1-ε*t)*(-0.7v3 + 0.054674v2*v3 + 0.00095888v2*v2*v3 - 0.064752v2*dv3 - 0.0007376v2*v2*dv3 - dv3) + ε*t*(-0.7v3 + 0.057618v2*v3 + 0.00068548v2*v2*v3 - 0.060546v2*dv3 - 0.00052729v2*v2*dv3 - dv3)
end

t_ = []
sol_ = [collect(traj0[1, vrbl[2]])[:, :]']
@showprogress for t0 in -200000:1000:300000
    prob = DE.DAEProblem(affined, collect(traj0[1, vrbl[1]]), last(eachrow(last(sol_)))[:], (t0, t0+1000), differential_vars = ones(Bool, length(vrbl[1])))
    sol = DE.solve(prob, Sundials.IDA(); initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-9)
    push!(t_, sol.t)
    push!(sol_, sol[1:3, :]')
end
t = vcat(t_...)
sol = vcat(sol_...)

points = arglmax(sol[:, 1])
scatter(t[points], sol[:, 1][points], ms = 1, msw = 0, color = :blue, alpha = 0.5, legend = :none, xformatter = t -> 1.7 + t/1000000, xlims = [-200000, 300000], ylims = [25, 40])
png("temp2")

