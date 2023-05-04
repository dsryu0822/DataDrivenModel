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
m, n = size(X)

@time using NoiseRobustDifferentiation
tvDX = tvdiff(X[1,:], 100, 1000, dx = 1)
plot(
    plot(X[1,:], title = "Data(WTI)"),
    plot(tvDX, title = "TVdiff"),
    plot(DX[1,:], title = "FDM"), layout = (3,1), size = (1000, 1000), legend = :none
)

ddprob = ContinuousDataDrivenProblem(X, DX)

@variables u[1:(size(X)[1])]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)
# basis = Basis(fourier_basis(u, 20), u)
# basis = Basis(sin_basis([i - j for j in u for i in u], 1), u)
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


using Interpolations

y = vec(rand(10))
cubic_f = cubic_spline_interpolation(1:10, y)
plot(1:0.1:10, cubic_f.(1:0.1:10), size = (800, 800)); scatter!(1:10, y)

Y_ = []
for i in 1:m
    y = vec(X[i, :])
    cubic_f = cubic_spline_interpolation(1:4610, y)
    push!(Y_, cubic_f.(1:0.1:4610))
end
Y = hcat(Y_...)'
DY, Y = fdiff(Y, stencil = 20, dt = 0.1, method = central_fdm)
ddprob = ContinuousDataDrivenProblem(Y, DY)

plot(1:0.1:4610, cubic_f.(1:0.1:4610), size = (800, 800))
scatter!(1:4610, y, msw = 0, ma = 0.5, ms = 2)