@time using OrdinaryDiffEq
@time using DataDrivenDiffEq
@time using DataDrivenSparse
@time using CSV, DataFrames

Oₘ = CSV.read("G:/DDM/atopy.txt", DataFrame, delim = '\t', header = ["P","B","D","R","K","G","kappa","t"])
dt = Oₘ.t[2] - Oₘ.t[1]

X = Matrix(Oₘ[1:1000:(end-1),1:3])
DX = ((Matrix(Oₘ[2:1000:end,1:3]) - Matrix(Oₘ[1:1000:(end-1),1:3])) / dt)'
# X = L[1:(end-1),:]
# DX = diff(L, dims = 1)'/0.001

ddprob = ContinuousDataDrivenProblem(X', DX)

@variables u[1:size(X')[1]]
u = DataDrivenDiffEq.scalarize(u)
basis1 = Basis(polynomial_basis(u, 4), u)
basis2 = Basis(fourier_basis(u, 400), u)
opt = STLSQ(0.1, 0)

ddsol1 = solve(ddprob, basis1, opt, options = DataDrivenCommonOptions(digits = 6))
@time ddsol2 = solve(ddprob, basis2, opt, options = DataDrivenCommonOptions(digits = 6))
# Basis(fourier_basis(u, 400), u) 54.902570 seconds - rss: 363

rss(ddsol1) / length(DX)
rss(ddsol2) / length(DX)

soleq1 = get_basis(ddsol1); get_parameter_values(soleq1)
soleq2 = get_basis(ddsol2); get_parameter_values(soleq2)

f̂1 = dynamics(soleq1);
P1 = get_parameter_values(soleq1)
f̂2 = dynamics(soleq2);
P2 = get_parameter_values(soleq2)

u0 = vec(X[1,:])
tspan = (0,10)
DDMsol1 = solve(ODEProblem(f̂1, u0, tspan, P1), RK4(), saveat = dt);
DDMsol2 = solve(ODEProblem(f̂2, u0, tspan, P2), RK4(), saveat = dt);
plot(
    plot(DDMsol2[2,:], ylabel = "B"),
    plot(DDMsol2(DDMsol2.t,Val{1})[2,:], ylabel = "dB/dt"),
    plot_title = "Result for Fourier Basis", layout = (2,1), legend = :none, lc = :black
)