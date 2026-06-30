include("../core/header.jl")
include("../core/factorio.jl")


trajargs = (; ticks = [], color = :black)
rcvdargs = (; ticks = [], color = :blue)
failargs = (; ticks = [], color = :red)
bfcnargs = (; yticks = [], msw = 0, ma = .5, ms = .5)

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    algae-zooplankton

''''''''''''''''''''''''''''''''''''''''''''''''''"""
vrbl = reverse(half(names(traj0)[Not(1)]))
cnfg = cookPI(vrbl; poly = 0:5)
traj0 = factory_algaezooplankton(DataFrame, 0.05, saveat = 0:1e-2:500)
traj1 = factory_algaezooplankton(DataFrame, 0.10, saveat = 0:1e-2:500)
plot(
    plot(traj0.A, traj0.Z),
    plot(traj1.A, traj1.Z),
)
f0 = SINDyPI(traj0, vrbl, cnfg; λ = 1e-8); f0 |> print
f1 = SINDyPI(traj1, vrbl, cnfg; λ = 1e-8); f1 |> print

rcvd0 = solve(DAEProblem(define(Function, f0), collect(traj0[1, vrbl[1]]), collect(traj0[1, vrbl[2]]), (0, 500), differential_vars = ones(Bool, length(vrbl[1]))), Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6) |> stack;
rcvd1 = solve(DAEProblem(define(Function, f1), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), (0, 500), differential_vars = ones(Bool, length(vrbl[1]))), Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6) |> stack;
plot(
    plot(eachcol(traj0[:, last(vrbl)])...; trajargs...),
    plot(eachcol(traj1[:, last(vrbl)])...; trajargs...),
    plot(eachrow(rcvd0)...; rcvdargs...),
    plot(eachrow(rcvd1)...; rcvdargs...),
    size = [400, 400]
)
define(f0, sigdigits = 5) |> print
define(f1, sigdigits = 5) |> print

function α_(t)
    maxα = 2.3
    return (maxα-abs(maxα - 5e-3t))
end
function affined(out, du, u, p, t)
    A, Z = u; dA, dZ = du
    α = α_(t)
    out11 = 0.6 + 1.4057A + 0.060348Z + 0.59668A*A - 0.52031A*Z + 0.0052994Z*Z - 0.13311A*A*A + 0.065296A*A*Z - 0.054688A*Z*Z - 0.0012047A*A*A*A - 0.0087065A*A*A*Z - 0.0010043A*A*Z*Z - 0.0058882A*Z*Z*Z - 2.5831e-5A*A*A*A*A - 0.0012826A*A*A*A*Z - 0.0012267A*A*A*Z*Z - 1.6761A*dA - 0.10058Z*dA - 0.015945A*A*dA - 0.17687A*Z*dA - 0.0088323Z*Z*dA - 0.00030997A*A*A*dA - 0.015391A*A*Z*dA - 0.014721A*Z*Z*dA - dA
    out21 = -0.15Z + 0.13274A*Z - 0.20172Z*Z + 0.017261A*A*Z - 0.35463A*Z*Z - 0.60229Z*Z*Z - 0.038357A*A*Z*Z + 0.52714A*Z*Z*Z - 0.006868Z*Z*Z*Z + 0.069043A*A*Z*Z*Z + 0.006868A*Z*Z*Z*Z - 1.7817A*dZ - 0.011447Z*dZ - 0.19178A*A*dZ - 0.019078A*Z*dZ - 4.0Z*Z*dZ - 7.127A*Z*Z*dZ - 0.045787Z*Z*Z*dZ - 0.76714A*A*Z*Z*dZ - 0.076312A*Z*Z*Z*dZ - dZ
    out12 = 0.6 + 1.402A + 0.062292Z + 0.58884A*A - 0.51625A*Z + 0.0060009Z*Z - 0.13499A*A*A + 0.070174A*A*Z - 0.055211A*Z*Z + 0.00035419A*A*A*A - 0.010423A*A*A*Z + 0.00019944A*A*Z*Z - 0.0066677A*Z*Z*Z - 0.00019449A*A*A*A*A - 0.0011739A*A*A*A*Z - 0.0013891A*A*A*Z*Z - 1.67A*dA - 0.10382Z*dA - 0.0069521A*A*dA - 0.18149A*Z*dA - 0.010002Z*Z*dA - 0.0023338A*A*A*dA - 0.014087A*A*Z*dA - 0.016669A*Z*Z*dA - dA
    out22 = -0.15Z + 0.15A*Z - 0.4Z*Z - 0.66667A*Z*Z - 0.6Z*Z*Z + 0.6A*Z*Z*Z - 1.6667A*dZ - 4.0Z*Z*dZ - 6.6667A*Z*Z*dZ - dZ

    out[1] = (1-α)*out11 + (α)*out12
    out[2] = (1-α)*out21 + (α)*out22
end
prob = DAEProblem(affined, collect(traj0[1, vrbl[1]]), [3,4], (0, 1000), differential_vars = ones(Bool, length(vrbl[1])))
sol = solve(prob, Sundials.IDA(); initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-36)
plot(
    # plot(sol.t, sol[1, :]),
    plot(sol.t, sol[2, :], ylims = [0, 5]),
    # plot(sol[1, :], sol[2, :]),
    plot([α_(t) for t in sol.t], sol[2, :], ylims = [0, 5]),
    plot(sol.t, [α_(t) for t in sol.t], ylims = [-.5, 2.5]),
)

result = DataFrame(t = sol.t, Z = sol[2, :], α = [α_(t) for t in sol.t], ud = [:blue, :red][[0; diff([α_(t) for t in sol.t]) .< 0] .+ 1])

plt_1 = plot(result.t, result.Z, color = result.ud, lw = 2, xticks = [0, 1000], yticks = [0, 4])
scatter!(plt_1, result.t[[1, end]], result.Z[[1, end]], color = result.ud[[1, end]], ms = 5, msc = :white, shape = [:circ, :rect])
plt_2 = plot(result.α, result.Z, color = result.ud, lw = 2, xticks = [0, 2.2], yticks = [0, 4])
scatter!(plt_2, result.α[[1, end]], result.Z[[1, end]], color = result.ud[[1, end]], ms = 5, msc = :white, shape = [:circ, :rect])
plt_3 = plot(result.t, result.α, color = result.ud, lw = 2, xticks = [0, 1000], yticks = [0, 2.2])

plot(plt_1, plt_2, plt_3)

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    pollinator-plant

''''''''''''''''''''''''''''''''''''''''''''''''''"""
traj0 = factory_pollinator(DataFrame, 0.74, ic = [.5, .5], saveat = 0:1e-2:100)
traj1 = factory_pollinator(DataFrame, 0.75, ic = [.5, .5], saveat = 0:1e-2:100)
plot(
    plot(traj0.P, traj0.A, color = :black),
    plot(traj1.P, traj1.A, color = :black),
    legend = :none, layout = (1, 2), size = [800, 400], lims = [0, 1.3]
)
vrbl = reverse(half(names(traj0)[Not(1)]))
cnfg = cookPI(vrbl; poly = 0:3)

f0 = SINDyPI(traj0, vrbl, cnfg; λ = 1e-5); f0 |> print
f1 = SINDyPI(traj1, vrbl, cnfg; λ = 1e-5); f1 |> print
define(f0, sigdigits = 5) |> print
define(f1, sigdigits = 5) |> print

function polli(out, du, u, p, t)
    P, A = u; dP, dA = du
    # ε = 1e-1
    α = (250-abs(250-t))/10
    out[1] = (1-α)*(0.3P - 1.0P*P + 2.1616P*A - 0.772P*P*A - 0.772A*dP - dP) +
               (α)*(0.3P - 1.0P*P + 2.1616P*A - 0.772P*P*A - 0.772A*dP - dP)  
    out[2] = (1-α)*(-0.44A + 1.4585P*A - 1.0A*A - 0.708P*A*A - 0.708P*dA - dA) + 
               (α)*(-0.45A + 1.4514P*A - 1.0A*A - 0.708P*A*A - 0.708P*dA - dA)
end

prob = DAEProblem(polli, collect(traj0[1, vrbl[1]]), collect(traj0[1, vrbl[2]]), (0, 1000), differential_vars = ones(Bool, length(vrbl[1])))
sol = solve(prob, Sundials.IDA(); initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-9)
plot(
    plot(sol.t, sol[1, :]),
    plot(sol.t, sol[2, :]),
    plot(sol[1, :], sol[2, :]),
    ylims = [0, Inf]
)

plot([(250-abs(250-t))/10 for t in 1:1000])

