using LinearAlgebra
using CSV, DataFrames
using Plots
include("../src/Utils.jl")

function hankelize(x, p, q)
    m = p*q
    H = zeros(q, p)
    for i ∈ 1:q
        for j ∈ 1:p
            H[i,j] = x[i+j-1]
        end
    end
    return H
end

dt = 0.001
r = 15

X = (CSV.read("test/lorenz.csv", DataFrame) |> Matrix)'[:, 1:(end-1)]

# X1 = X[:, 1:(end-1)]
# X2 = circshift(X, (0,-1))[:, 1:(end-1)]

# U, Σ, V = svd(X1); Σ = diagm(Σ)

# Uᵣ = U[:, 1:r]
# Σᵣ = Σ[1:r, 1:r]
# Vᵣ = V[:, 1:r]
# Ã = Uᵣ' * X2 * Vᵣ / Σᵣ
# D, Wᵣ = eigen(Ã); D = diagm(D)
# Φ = X2 * Vᵣ / Σᵣ * Wᵣ

# reconstructed = real.(Φ * X1)
# plot(X1[1,(50000:60000)], label = "Simulation")
# plot!(-reconstructed[1,(50000:60000)], label = "Reconstructed")
# plot(eachrow(X1)...),
# plot(eachrow(-reconstructed)...)

x = X[1,:]
H = hankelize(x, length(x) - 100, 100)

(U, Σ, V) = svd(H) # U * diagm(Σ) * V' ≈ H
r = min(r, size(U, 2)); V = V[:,1:r]'
# plot(V[:,1], V[:,2], V[:,3])
# plot(V[1:25000,1], size = (1000, 300))
# plot(x[1:100000], size = (1000, 300))

dV, V = fdiff(V, stencil = 4, dt = dt, method = central_fdm); dV = dV[:,1:end-1]; V = V[:,1:end-1]
# normΘ = norm.(eachrow(V))
# V = V ./ normΘ
plot(
    plot(x[1:25000], label = "x"),
    plot(Σ[  1] * V[1,1:25000], label = "v1"),
    plot(Σ[  2] * V[2,1:25000], label = "v2"),
    plot(Σ[  3] * V[3,1:25000], label = "v2"),
    plot(Σ[  4] * V[4,1:25000], label = "v2"),
    plot(Σ[  5] * V[5,1:25000], label = "v2"),
    plot(Σ[  6] * V[6,1:25000], label = "v2"),
    plot(Σ[  7] * V[7,1:25000], label = "v2"),
    plot(Σ[  8] * V[8,1:25000], label = "v2"),
    plot(Σ[end] * V[end,1:25000], label = "vr"),
    layout = (:,2)
)
plot(
    # plot(x),
    plot((x[1:199900] .- reduce(+, Σ[1:r] .* eachrow(V)))[1:25000])
)

@time using DataDrivenDiffEq
@time using DataDrivenSparse
ddprob = ContinuousDataDrivenProblem(V, dV)

@variables ξ[1:r]
ξ = DataDrivenDiffEq.scalarize(ξ)
basis = Basis(polynomial_basis(ξ, 1), ξ)

opt = STLSQ(1.1)
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 6))
soleq = get_basis(ddsol);
print(soleq, true)
get_parameter_map(soleq)

f̂ = dynamics(soleq)
P = get_parameter_values(soleq)
# f̂(u0, P, 0)

# @time using OrdinaryDiffEq
# DDM = solve(ODEProblem(f̂, V[:,1], (0,25), P), RK4(), saveat = dt)
# plot(
#     plot(DDM, idxs = (0,1), color = :blue),
#     plot(DDM, idxs = (1,2,3), color = :blue),
#     plot(DDM, idxs = (0,15), color = :blue),
#     size = (600, 600), layout = (3,1)
# )
# V