include("dataload.jl")

@time using OrdinaryDiffEq
@time using DataDrivenSparse

selected = [:WTI유, :구리, :금]
# X = Matrix(trng[:, selected])'
X = Matrix(trng[:, Not(:t)])'
X = Array(sol)
# DX = diff(X, dims = 2); X = X[:, 1:(end-1)]; DX
# DX = sol(sol.t,Val{1})[1:3,:]
DX = fdiff(X; order = 2, dt = dt); X = X[:, 2:(end-0)]; DX
# diff(X, dims = 2) / dt

# ddprob = DataDrivenProblem(sol)
p = sol.prob.p
t = sol.t
DX = similar(X)
if DiffEqBase.isinplace(sol.prob.f)
    foreach(enumerate(eachcol(X))) do (i, xi)
        @views sol.prob.f(DX[:, i], xi, p, t[i])
    end
else
    foreach(enumerate(eachcol(X))) do (i, xi)
        DX[:, i] .= sol.prob.f(xi, p, t[i])
    end
end
DX

ddprob = ContinuousDataDrivenProblem(X, DX)

@variables u[1:size(X)[1]]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)

opt = STLSQ(exp10.(-3:0.1:-1), 10)
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 5))
soleq = get_basis(ddsol);
get_parameter_map(soleq)
print(soleq, true) # soleq[1].rhs
is_converged(ddsol)

f̂ = dynamics(soleq);
P = get_parameter_values(soleq)

u0 = collect(trng[end, Not(1)])
tend = 250
tspan = (30, tend)
dt = 0.01
DDM = solve(ODEProblem(f̂, u0, tspan, P), RK4(), saveat = dt);

@time using Plots; default(fontfamily = "맑은 고딕")
plt_ = []
for k in 1:size(X)[1]
    ptemp = plot(DDM, idxs = (0,k), legend = :none)
    k = k+1
    plot!(ptemp, 0:30, trng[(end-30):end, k], xlims = (0,tend), color = :black, title = names(data)[k])
    plot!(ptemp, 30:tend, test[1:(tend-29), k], color = :black)
    push!(plt_, ptemp)
end
plot(plt_..., size = (1000,1000), layout = (:,2))

DDM = solve(ODEProblem(f̂, u0, (0,100), P), RK4(), saveat = dt)

plot(
    plot(sol, vars = (1, 2, 3), color = :black),
    plot(DDM, vars = (1, 2, 3), color = :blue),
    size = (600, 600)
)
