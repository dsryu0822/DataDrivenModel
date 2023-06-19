using CSV, DataFrames, LinearAlgebra, Clustering, Plots

function RK4(f::Function,v::AbstractVector, h=10^(-2))
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4)
end


m = 10^(-3)
μ = 10^(-6)
R = 22
C = 47μ
L = 20m
T = 400μ
γ = 11.75238
η = 1309.524
E = 53.500001

Vr(t) = γ + η * (mod(t, T))
function buck(v)
    V, I, t = v

    V̇ = - V/(R*C) + I/C
    İ = - (V/L) + ifelse(V < Vr(t), E/L, 0)
    return [V̇, İ, 1]
end

u0 = [12.3, 0.55, 0.0]
u_ = [u0]
∇_ = []
dt = 0.00001; tend = 0.25
for t in dt:dt:tend
    push!(∇_, buck(u_[end]))
    push!(u_, RK4(buck, u_[end], dt))
end
push!(∇_, buck(u_[end]))
U = stack(u_)[Not(end), :]
∇ = stack(∇_)[Not(end), :]
DATA = DataFrame(
    [collect(0:dt:tend)'; U; ∇; Vr.(0:dt:tend)']'
    , ["t", "V", "I", "dV", "dI", "Vr"])

Y = select(DATA, [:dV, :dI]) |> Matrix .|> Float32
X = select(DATA, [ :V,  :I]) |> Matrix .|> Float32
XY = [X Y]

## Clustering
dbs = dbscan(col_normalize(Y)', 0.01); nsubsys = length(dbs.clusters); println(nsubsys, " clusters found!")
plot(DATA.V, DATA.I, color = dbs.assignments, alpha = 0.5)
bit_ = [dbs.assignments .== k for k in 1:nsubsys]
Θ_ = [poly_basis(X[s,:], 2)     for s in bit_]
Y_ = [Y[vcat(findall(s)...),:] for s in bit_]
Ξ_ = [STLSQ(Θ_[s], Y_[s], 0.01) for s in 1:nsubsys]

[
    0        0
    -1/(R*C) -1/L
    1/C      0
    0        0
    0        0
    0        0
]
Ξ_[1]


[
    0        E/L
    -1/(R*C) -1/L
    1/C      0
    0        0
    0        0
    0        0
]
Ξ_[2]
