include("../core/header.jl")
include("tippingutils.jl")

using DifferentialEquations
import Sundials
import DiffEqBase

function factory_pollinator(κ::Number; ic = [0., 1, 1], tspan = 0:1e-2:10)
      α, β,  _h,   γp,   γA = (
    0.3, 1, 0.4, 1.93, 1.77 )
    function sys(v::AbstractVector)
        _,P,A = v

        dt = 1
        dP = P*(α - β*P + frac(γp*A, 1 + _h*γp*A))
        dA = A*(α - κ - β*A + frac(γA*P, 1 + _h*γA*P))
        return [dt, dP, dA]
    end
        
    len_t_ = length(tspan); h = tspan.step.hi
    v = ic; DIM = length(v)
    traj = zeros(2+len_t_, 2DIM); traj[1, 1:DIM] = v

    tk = 0
    while tk ≤ len_t_
        t,_,_ = v
        v, dv = RK4(sys, v, h)

        if t ≥ first(tspan)
            tk += 1
            traj[tk+1,         1:DIM ] =  v
            traj[tk  , DIM .+ (1:DIM)] = dv
        end
    end
    return traj[1:(end-2), :]
end
factory_pollinator(T::Type, args...; kargs...) =
DataFrame(factory_pollinator(args...; kargs...), ["t", "P", "A", "dt", "dP", "dA"])[:, Not(:dt)]

traj0 = factory_pollinator(DataFrame, 0.74, ic = [0, .5, .5], tspan = 0:1e-2:100)
traj1 = factory_pollinator(DataFrame, 0.75, ic = [0, .5, .5], tspan = 0:1e-2:100)
traj2 = factory_pollinator(DataFrame, 0.90, ic = [0, .5, .5], tspan = 0:1e-2:100)
plot(
    plot(traj0.P, traj0.A, color = :black),
    plot(traj1.P, traj1.A, color = :black),
    plot(traj2.P, traj2.A, color = :black),
    legend = :none, layout = (1, 3), size = [1200, 400], lims = [0, 1.3]
)

vrbl = reverse(half(names(traj0)[Not(1)]))
cnfg = cookPI(vrbl; poly = 0:3)

f0 = SINDyPI(traj0, vrbl, cnfg; λ = 1e-5); f0 |> print
f1 = SINDyPI(traj1, vrbl, cnfg; λ = 1e-5); f1 |> print
f2 = SINDyPI(traj2, vrbl, cnfg; λ = 1e-5); f2 |> print

define(f2) |> print

tspan = (0, 100)
sargs = (; reltol = 1e-6, initializealg = DiffEqBase.BrownFullBasicInit())
prob0 = DAEProblem(define(Function, f0), collect(traj0[1, vrbl[1]]), collect(traj0[1, vrbl[2]]), tspan, differential_vars = ones(Bool, length(vrbl[1])))
sol0 = solve(prob0, Sundials.IDA(); sargs...)

prob1 = DAEProblem(define(Function, f1), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), tspan, differential_vars = ones(Bool, length(vrbl[1])))
sol1 = solve(prob1, Sundials.IDA(); sargs...)

prob2 = DAEProblem(define(Function, f2), collect(traj2[1, vrbl[1]]), collect(traj2[1, vrbl[2]]), tspan , differential_vars = ones(Bool, length(vrbl[1])))
sol2 = solve(prob2, Sundials.IDA(); sargs...)

plot(
    plot(eachrow(stack(sol0.u))..., color = :blue),
    plot(eachrow(stack(sol1.u))..., color = :blue),
    plot(eachrow(stack(sol2.u))..., color = :blue),
    legend = :none, layout = (1, 3), size = [1200, 400], lims = [0, 1.3]
)

α_ = 0:.1:26
vrtl1 = []
vrtl2 = []
@showprogress for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg[1])
    probg = DAEProblem(define(Function, g), collect(traj0[1, vrbl[1]]), collect(traj0[1, vrbl[2]]), (0, 1000), differential_vars = ones(Bool, length(vrbl[1])))
    solg = solve(probg, Sundials.IDA(); sargs...)
    traj = stack(solg.u)'[end, :]
    push!(vrtl1, traj[1])
    push!(vrtl2, traj[2])
end
plot(
    plot(α_, vrtl1, ms = 1, msw = 0, color = :black, alpha = 0.5, legend = :none, ylims = [-0, 2], ylabel = "P"),
    plot(α_, vrtl2, ms = 1, msw = 0, color = :black, alpha = 0.5, legend = :none, ylims = [-1, 2], ylabel = "A"),
    layout = (2, 1), size = [400, 400]
)
png("temp")

