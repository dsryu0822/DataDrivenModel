include("dataload.jl")
include("Utils.jl")

@time using OrdinaryDiffEq
@time using DataDrivenDiffEq
@time using DataDrivenSparse
@time using Plots; default(fontfamily = "")

selected = [:WTI유, :구리]
X = Matrix(trng[:, Not(:t)])'
# X = Matrix(trng[:, selected])'
# DX = diff(X, dims = 2); X = X[:, 1:(end-1)]; DX
DX, X = fdiff(X, stencil = 2, dt = 1, method = central_fdm)

ddprob = ContinuousDataDrivenProblem(X, DX)

@variables u[1:size(X)[1]]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)
# basis = Basis(fourier_basis(u, 20), u)
# basis = Basis(cos_basis([i - j for j in u for i in u], 1), u)
# basis = Basis([fourier_basis(u, 20); polynomial_basis(u, 2)], u)

opt = STLSQ(10^(-6), 10)
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 6))
soleq = get_basis(ddsol);
get_parameter_map(soleq)
print(soleq, true) # soleq[1].rhs
is_converged(ddsol)

f̂ = dynamics(soleq);
P = get_parameter_values(soleq)

u0 = collect(trng[end, Not(1)])
tend = 250; tspan = (30, tend)
dt = 1
DDM = solve(ODEProblem(f̂, u0, tspan, P), RK4(), saveat = dt)

plt_ = []
for k in 1:size(X)[1]
    ptemp = plot(DDM, idxs = (0,k), legend = :none)
    k = k+1
    plot!(ptemp, 0:30, trng[(end-30):end, k], xlims = (0,tend), color = :black, title = names(data)[k])
    plot!(ptemp, 30:tend, test[1:(tend-29), k], color = :black)
    push!(plt_, ptemp)
end
plot(plt_..., size = (1000,1000), layout = (:,2))

plot(vec(sum(X, dims = 1)))
plot(vec(sum(DX, dims = 1)))

histogram(vec(sum(DX, dims = 1)))