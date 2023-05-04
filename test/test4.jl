@time using OrdinaryDiffEq
@time using DataDrivenSparse
@time using DataDrivenDiffEq
@time using Plots

function f(vec, p, t)
    x, y, z = vec

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z

    return [ẋ, ẏ, ż]
end
p,t = Nothing, Nothing
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
dt = 0.01
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, RK4(), saveat = dt)

X = Array(sol)
DX1 = hcat(sol.prob.f.(eachcol(X), p, t)...)
DX2 = Array(sol(sol.t, Val{1}))

ddprob1 = ContinuousDataDrivenProblem(X, DX1)
ddprob2 = ContinuousDataDrivenProblem(X, DX2)

@variables u[1:size(X)[1]]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)
opt = STLSQ(exp10.(-5:0.1:-1), 0)

ddsol1 = solve(ddprob1, basis, opt, options = DataDrivenCommonOptions(digits = 6))
ddsol2 = solve(ddprob2, basis, opt, options = DataDrivenCommonOptions(digits = 6))

soleq1 = get_basis(ddsol1); get_parameter_values(soleq1)
soleq2 = get_basis(ddsol2); get_parameter_values(soleq2)

f̂1 = dynamics(soleq1);
P1 = get_parameter_values(soleq1)
f̂2 = dynamics(soleq2);
P2 = get_parameter_values(soleq2)

DDMsol1 = solve(ODEProblem(f̂1, u0, tspan, P1), RK4(), saveat = dt);
DDMsol2 = solve(ODEProblem(f̂2, u0, tspan, P2), RK4(), saveat = dt);

plot(
    plot(sol, idxs = (1,2,3), title = "Ground Truth"),
    plot(DDMsol1, idxs = (1,2,3), title = "Exact"),
    plot(DDMsol2, idxs = (1,2,3), title = "Order 1"),
    size = (600, 600)
)

rss(ddsol1)
rss(ddsol2)

using FiniteDifferences

central_fdm(5, 1)
backward_fdm(2, 1)

include("../src/Utils.jl")
fdiff(X, stencil = 2, dt = dt)
DX2[:, 2:end]

fdiff(X, stencil = 5, dt = dt)
DX1[:, 5:end]
DX2[:, 5:end]

sum(abs, DX1[:, 5:end] - fdiff(X, order = 5, dt = dt))
sum(abs, DX1[:, 5:end] - DX2[:, 5:end])

# https://symbolics.juliasymbolics.org/dev/manual/sparsity_detection/
bart = Dict(get_parameter_map(soleq1)) |> keys
foo = Symbolics.arguments(Symbolics.arguments(soleq1[2].rhs)[3])[1]
foo in bart