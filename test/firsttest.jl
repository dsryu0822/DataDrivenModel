@time using DataDrivenDiffEq
@time using ModelingToolkit
@time using OrdinaryDiffEq
@time using DataDrivenSparse
@time using LinearAlgebra
@time using Plots

# Create a test problem
function lorenz(u, p, t)
    x, y, z = u

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z
    return [ẋ, ẏ, ż]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
dt = 0.1
prob = ODEProblem(lorenz, u0, tspan)
# sol = solve(prob, Tsit5(), dense = true) # sol(sol.t,Val{4}); sol(sol.t,Val{5})
sol = solve(prob, RK4(), saveat = dt) # Tsit5() is more efficient https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Explicit-Runge-Kutta-Methods
plot(
    plot(sol, idxs = (0, 1), color = :black),
    plot(sol, idxs = (0, 2), color = :black),
    plot(sol, idxs = (0, 3), color = :black),
    plot(sol, idxs = (1, 2, 3), color = :black)
)
# How get the derivative from ODEsolution
# https://github.com/SciML/DataDrivenDiffEq.jl/blob/8acc287468ef13be134ebd2d85d273bc0f456cfb/src/problem/type.jl#L486
@assert sol(sol.t,Val{1})[1] == ((sol(dt) - sol(0)) / dt)

## Start the automatic discovery
ddprob = DataDrivenProblem(sol)

@variables x y z
u = [x; y; z]
basis = Basis(polynomial_basis(u, 2), u)
# @variables t x(t) y(t) z(t)
# u = [x; y; z]
# basis = Basis(polynomial_basis(u, 2), u, iv = t)

opt = STLSQ(exp10.(-5:0.1:-1))
# https://github.com/SciML/DataDrivenDiffEq.jl/blob/cae0c79e0d7b9eef48f9350e01f69b0798b447bb/lib/DataDrivenSparse/src/algorithms/STLSQ.jl#L97-L120
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 4))
soleq = get_basis(ddsol)
println(soleq) # soleq[1].rhs
get_parameter_map(soleq) # Equlvalent to get_parameter_values(soleq)
ddsol.residuals # Equlvalent to rss(ddsol)