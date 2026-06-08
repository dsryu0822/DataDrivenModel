include("../core/header.jl")
import DifferentialEquations as DE

function thomas(du, u, p, t)
    x, y, z = u
    b = p[1]

    du[1] = sin(y) - b*x
    du[2] = sin(z) - b*y
    du[3] = sin(x) - b*z
    return du
end

bfcn = Dict()
@showprogress @threads for b in .10:1e-4:.24
    sol = DE.solve(DE.ODEProblem(thomas, rand(3), (0, 2000), [b]), DE.RK4(), dt = 1e-2, adaptive=false)
    x_ = sol[1, sol.t .> 1000]
    bfcn[b] = x_[arglmax(x_)]
end
scatter([[fill(k, length(v)) for (k,v) in bfcn]...;], [values(bfcn)...;], msw = 0, ms = .5, color = :black, ma = .5, legend = :none, xflip = true, xticks = .1:.01:.24, ygrid = false, xlabel = L"b", ylims = [-10, 10], xlims = [.1, .24])
png("G:/bifurcation.png")


function factory_thomas(b::Number; ic = rand(3), tspan = 0:1e-2:2000)
    function thomas(du, u, p, t)
        x, y, z = u
        b = p[1]

        du[1] = sin(y) - b*x
        du[2] = sin(z) - b*y
        du[3] = sin(x) - b*z
        return du
    end
    t0 = first(tspan)
    sol = DE.solve(DE.ODEProblem(thomas, ic, (0, last(tspan)), [b]), DE.RK4(), dt = tspan.step.hi, adaptive=false)
    matrix = Matrix([sol.t'; sol[:, :]; stack([thomas(zeros(3), u, [b], 0) for u in sol.u])]')
    return matrix[sol.t .≥ t0, :][1:end-1, :]
end
factory_thomas(T::Type, args...; kargs...) =
DataFrame(factory_thomas(args...; kargs...), ["t", "x", "y", "z", "dx", "dy", "dz"])

traj0 = factory_thomas(DataFrame, .140, tspan = 1000:1e-2:2000)
plt_pp_1 = plot(traj0.x, traj0.y, traj0.z, color = :black)
traj1 = factory_thomas(DataFrame, .150, tspan = 1000:1e-2:2000)
plt_pp_2 = plot(traj1.x, traj1.y, traj1.z, color = :black)

vrbl = reverse(half(names(traj0)[2:end]))
# cnfg = cook(last(vrbl); poly = 0:5)
cnfg = cook(last(vrbl); poly = 0:1, trig = [1], format = sin)

@time f0 = SINDy(traj0, vrbl, cnfg; λ = 1e-3); f0 |> print
@time f1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); f1 |> print

prob0 = DE.ODEProblem(define(Function, f0), collect(traj0[1, vrbl[2]]), (0, 2000))
sol0 = DE.solve(prob0, DE.RK4(), dt = 1e-2, adaptive=false);
plt_pp_3 = plot(eachrow(stack(sol0.u))...; color = :blue)
prob1 = DE.ODEProblem(define(Function, f1), collect(traj1[1, vrbl[2]]), (0, 2000))
sol1 = DE.solve(prob1, DE.RK4(), dt = 1e-2, adaptive=false);
plt_pp_4 = plot(eachrow(stack(sol1.u))...; color = :blue)
plot(plt_pp_1, plt_pp_2, plt_pp_3, plt_pp_4, legend = :none, size = [600, 600])

α_ = -4:1e-2:10
f_ = []
for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg, method = "SINDy")

    push!(f_, define(Function, g; fname = "f_$k", sigdigits = 20))
end
define(f_[1]) |> println

bfcn = Dict()
@showprogress for k in eachindex(α_)
    try
        sol = DE.solve(DE.ODEProblem(f_[k], rand(3), (0, 2000), []), DE.RK4(), dt = 1e-2, adaptive=false)
        x_ = sol[1, sol.t .> 1000]
        bfcn[α_[k]] = x_[arglmax(x_)]
    catch e
        continue
    end
end
filter(x -> !isempty(second(x)), bfcn)
for (k, v) in bfcn
    if v |> isempty
        delete!(bfcn, k)
    end
end
scatter([[fill(k, length(v)) for (k,v) in bfcn]...;], [values(bfcn)...;], msw = 0, ms = .5, color = :blue, ma = .5, legend = :none, xflip = true, xlims = [-4, 10], xticks = -4:1:10, ylims = [-10, 10])
