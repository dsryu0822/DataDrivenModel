include("../core/DDM.jl")
include("../core/header.jl")
include("tippingutils.jl")

import DifferentialEquations as DE
import Sundials
import DiffEqBase

function lorenz!(du, u, p, t)
    σ, β = p
    x, y, z, ρ = u
    du[1] = σ*(y - x)
    du[2] = x*(ρ - z) - y
    du[3] = x*y - β*z
    du[4] = -1e-2
end
prob = DE.ODEProblem(lorenz!, [1, 1, 1, 230], (0, 3000), [10, 8/3])
sol = DE.solve(prob, DE.RK4(), saveat = 0:1e-3:3000)

bit_xm = arglmax(sol[3, :])
scatter(sol.t[bit_xm], sol[3, :][bit_xm], msw = 0, ylims = [220, 320], ms = 1, color = :black, legend = :none, alpha = 0.5, xformatter = t -> 230 - t/100, xlims = [0, 3000])

# hrzn = []
# vrtl = []
# for ρ in 230:(-0.1):200
#     traj = factory_lorenz(DataFrame, ρ, tspan = [0., 50.], dt = 1e-3)[40000:end, :]
#     points = arglmax(traj.z)
#     push!(vrtl, traj.z[points])
#     push!(hrzn, repeat([ρ], length(points)))
# end
# scatter([hrzn...;], [vrtl...;], xflip = true, xlims = [200, 230], ylims = [220, 320], msw = 0, color = :black, ms = 1, legend = :none, alpha = 0.5)
# png("3")

traj0 = factory_lorenz(DataFrame, 230, tspan = [0., 50.], dt = 1e-3)[40000:end, :]
traj1 = factory_lorenz(DataFrame, 225, tspan = [0., 50.], dt = 1e-3)[40000:end, :]
plot(
    plot(traj0.x, traj0.y, traj0.z, ticks = [], color = :black),
    plot(traj1.x, traj1.y, traj1.z, ticks = [], color = :black),
    legend = :none
)


    vrbl = reverse(half(names(traj0)))
    cnfg = cookPI(vrbl; poly = 0:2)

    @time f0 = SINDyPI(traj0, vrbl, cnfg, λ = 0.01); f0 |> print
    @time f1 = SINDyPI(traj1, vrbl, cnfg, λ = 0.01); f1 |> print

    tspan = (0, 50)
    sargs = (; reltol = 1e-6, initializealg = DiffEqBase.BrownFullBasicInit(), saveat = 0:1e-3:100)
    prob0 = DE.DAEProblem(define(Function, f0), collect(traj0[1, vrbl[1]]), collect(traj0[1, vrbl[2]]), tspan, differential_vars = ones(Bool, length(vrbl[1])))
    sol0 = DE.solve(prob0, Sundials.IDA(); sargs...)

    define(f1) |> print

    prob1 = DE.DAEProblem(define(Function, f1), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), tspan, differential_vars = ones(Bool, length(vrbl[1])))
    sol1 = DE.solve(prob1, Sundials.IDA(); sargs...)

    plot(
        plot(eachrow(stack(sol0.u))..., ticks = [], color = :blue),
        plot(eachrow(stack(sol1.u))..., ticks = [], color = :blue),
        legend = :none, layout = (1, 3), size = [1200, 400]
    )

α_ = 0:1e-1:6
hrzn = []
vrtl = []
@showprogress for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg[1])
    probg = DE.DAEProblem(define(Function, g), collect(traj0[1, vrbl[2]]), collect(traj0[1, vrbl[1]]), (0, 100), differential_vars = ones(Bool, length(vrbl[1])))
    solg = DE.solve(probg, Sundials.IDA(); sargs...)
    traj = solg[3, :][solg.t .> 90]
    points = arglmax(traj[:, 1])
    push!(vrtl, traj[points, 1])
    push!(hrzn, repeat([α], length(points)))
end
trajg = DataFrame(stack(solg.u)', [:R, :C, :P])
scatter(hrzn, vrtl, ms = 1, msw = 0, color = :blue, alpha = 0.5, legend = :none, xlims = [0, 6], ylims = [220, 320])
png("temp")

# function affinize(f0, f1)
#     V1 = replace(define(f0), "out" => "out0") * replace(define(f1), "out" => "out1")
#     V2 = replace(V1, r"end.+= du"s => "", "f(out0"=> "f(out", "end" => "out = (1-t)*out0 + t*out1\nend")
#     V4 = replace(V2, "u, p, t)" => "u, p, t)\nout0, out1 = deepcopy(out), deepcopy(out)\n")
#     return V4
# end
# affinize(f0, f1) |> print

# define(f0) |> print
# define(f1) |> print

function affined(out, du, u, p, t)
x, y, z = u; dx, dy, dz = du
ε = 1e-3
out[1] = (1 - ε*t)*(-10.0x + 10.0y)         + ε*t*(-10.0x + 10.0y) - dx
out[2] = (1 - ε*t)*(230.0x - 1.0y - 1.0x*z) + ε*t*(225.0x - 1.0y - 1.0x*z) - dy
out[3] = (1 - ε*t)*(-2.6667z + 1.0x*y)      + ε*t*(-2.6667z + 1.0x*y) - dz
end

t_ = []
sol_ = [collect(traj0[1, vrbl[2]])[:, :]']
@showprogress for t0 in 0:100:5900
    prob = DE.DAEProblem(affined, collect(traj0[1, vrbl[1]]), last(eachrow(last(sol_)))[:], (t0, t0+100), differential_vars = ones(Bool, length(vrbl[1])))
    sol = DE.solve(prob, Sundials.IDA(); initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6)
    push!(t_, sol.t)
    push!(sol_, sol[1:3, :]')
end
t = vcat(t_...)
sol = vcat(sol_...)

points = arglmax(sol[:, 3])
scatter(t[points], sol[:, 3][points], ms = 1, msw = 0, color = :blue, alpha = 0.5, legend = :none, ylims = [220, 320], xformatter = t -> 230 - t/200, xlims = [0, 6000])
png("temp2")


plot(t[1:100:end], sol[1:100:end, 3], xlabel = "t", ylabel = "z", legend = :none)

scatter(up(sol[:, 3][points]), dw(sol[:, 3][points]))