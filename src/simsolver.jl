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

p,t = Nothing, Nothing
u0 = [1; 0; 0]
u0 = [-8; 8; 27]
tspan = (0.0, 200.0)
dt = 0.001
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, RK4(), saveat = dt)
plot(sol, idxs = (1, 2, 3), color = :black)
@assert sol(sol.t,Val{1})[1] == ((sol(dt) - sol(0)) / dt)

ddprob = DataDrivenProblem(sol)

using FiniteDifferences
include("Utils.jl")
X = Array(sol)
DX, X = fdiff(X, stencil = 2, dt = dt, method = central_fdm)
# DX = hcat(sol.prob.f.(eachcol(X), p, t)...)
ddprob = ContinuousDataDrivenProblem(X, DX)

Y_ = []
for i in 1:3
    y = sol[i,:]
    cubic_f = cubic_spline_interpolation(0:dt:100, y)
    push!(Y_, cubic_f.(0:(0.0001):100))
end
Y = hcat(Y_...)'
DY, Y = fdiff(Y, stencil = 4, dt = 0.0001, method = central_fdm)
ddprob = ContinuousDataDrivenProblem(Y, DY)

xyz = @variables x y z
basis = Basis(polynomial_basis(xyz, 2), xyz)

opt = STLSQ(exp10.(-6:1:-1))
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 6))
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

# badpoints = []
for k in 9000:(-10):1
    residuals = [sum(abs2, f̂(X[:,j], P, 0) - DX[:,j]) for j in axes(X)[2]]
    # residuals[badpoints] .= 0
    # push!(badpoints, argmax(residuals))
    # taboo = Not(badpoints)
    taboo = ordinalrank(residuals) .< k
    newddprob = ContinuousDataDrivenProblem(X[:, taboo], DX[:, taboo])
    newddsol = solve(newddprob, basis, opt, options = DataDrivenCommonOptions(digits = 4))
    newsoleq = get_basis(newddsol)
    f̂ = dynamics(newsoleq)
    P = get_parameter_values(newsoleq)
        println(get_parameter_map(newsoleq))
        println("k = ", k, "\trss: ", rss(newddsol))
    # print(newsoleq, true)
end

using StatsBase
x = rand(10)

