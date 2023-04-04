@time using DataDrivenDiffEq
@time using ModelingToolkit
@time using OrdinaryDiffEq
@time using DataDrivenSparse
@time using LinearAlgebra
@time using Plots

# Create a test problem
function f(vec, p, t)
    x, y, z, u, v, w = vec

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z

    u̇ = - v - w
    v̇ = u + 0.1 * v
    ẇ = 0.1 + w * (u - 14)

    return [ẋ, ẏ, ż, u̇, v̇, ẇ]
end

u0 = [1.0; 0.0; 0.0; 1.0; 1.0; 1.0]
tspan = (0.0, 100.0)
dt = 0.01
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, RK4(), saveat = dt)
plot(
    plot(sol, idxs = (1, 2, 3), color = :black),
    plot(sol, idxs = (4, 5, 6), color = :black)
)
@assert sol(sol.t,Val{1})[1] == ((sol(dt) - sol(0)) / dt)

ddprob = DataDrivenProblem(sol)

vec = @variables x y z u v w
basis = Basis(polynomial_basis(vec, 2), vec)

opt = STLSQ(exp10.(-5:0.1:-1))
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 8))
soleq = get_basis(ddsol)
print(soleq, true) # soleq[1].rhs
get_parameter_map(soleq) # Equlvalent to get_parameter_values(soleq)
ddsol.residuals # Equlvalent to rss(ddsol)


f̂ = dynamics(soleq)
P = get_parameter_values(soleq)
# f̂(u0, P, 0)

DDM = solve(ODEProblem(f̂, u0, tspan, P), RK4(), saveat = dt)
plot(
    plot(sol, idxs = (1, 2, 3), color = :black),
    plot(sol, idxs = (4, 5, 6), color = :black, zlabel = "SimulationModel"),
    plot(DDM, idxs = (1, 2, 3), color = :blue),
    plot(DDM, idxs = (4, 5, 6), color = :blue, zlabel = "DataDrivenModel"),
    size = (600, 600)
)