# function f(vec, p, t)
#     x, y, z, u, v, w = vec

#     ẋ = 10.0 * (y - x)
#     ẏ = x * (28.0 - z) - y
#     ż = x * y - (8 / 3) * z

#     u̇ = - v - w
#     v̇ = u + 0.1 * v
#     ẇ = 0.1 + w * (u - 14)

#     return [ẋ, ẏ, ż, u̇, v̇, ẇ]
# end
function f(vec, p, t)
    x, y, z = vec

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z

    return [ẋ, ẏ, ż]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
dt = 0.01
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, RK4(), saveat = dt)
plot(sol, vars = (1, 2, 3), color = :black)
@assert sol(sol.t,Val{1})[1] == ((sol(dt) - sol(0)) / dt)

ddprob = DataDrivenProblem(sol)

vec = @variables x y z
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
    plot(sol, vars = (1, 2, 3), color = :black),
    plot(DDM, vars = (1, 2, 3), color = :blue),
    size = (600, 600)
)