include("dataload.jl")
include("Utils.jl")

@time using OrdinaryDiffEq
@time using DataDrivenDiffEq
@time using DataDrivenSparse
@time using Plots; default(plot_titlefontfamily = "굴림체")

selected = [:WTI유, :구리, :금]
# X = Matrix(trng[:, selected])'
X = Matrix(trng[:, Not(:t)])'
# DX = diff(X, dims = 2); X = X[:, 1:(end-1)]; DX
DX, X = fdiff(X, stencil = 4, dt = 1, method = backward_fdm)

ddprob = ContinuousDataDrivenProblem(X, DX)

@variables u[1:size(X)[1]]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)

opt = STLSQ(10^(-3), 10)
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 4))
soleq = get_basis(ddsol);
get_parameter_map(soleq)
print(soleq, true) # soleq[1].rhs
is_converged(ddsol)

f̂ = dynamics(soleq);
P = get_parameter_values(soleq)

u0 = collect(trng[end, Not(1)])
tend = 250
tspan = (30, tend)
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