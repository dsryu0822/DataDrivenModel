using CSV, DataFrames
using Plots
using Combinatorics

function libararize(X, K = 1)
    mtrx_X = Matrix(X)
    dim = size(X, 2)
    ansatz = []

    for k in 0:K
        for case = collect(multiexponents(dim,k))
            push!(ansatz, prod(mtrx_X .^ case', dims = 2))
        end
    end
    return hcat(ansatz...)
end

Oₘ = CSV.read("G:/DDM/atopy.txt", DataFrame, delim = '\t', header = ["P","B","D","R","K","G","kappa","t"]) #, limit = 10000)
dt = Oₘ.t[2] - Oₘ.t[1]
Oₘ.B .*= 100

# Oₘ.P ./= sqrt(sum(abs2,Oₘ.P))
# Oₘ.B ./= sqrt(sum(abs2,Oₘ.B))
# Oₘ.D ./= sqrt(sum(abs2,Oₘ.D))


# plot(
#     plot(Oₘ.t[50000:60000], Oₘ.P[50000:60000], title = "P"),
#     plot(Oₘ.t[50000:60000], Oₘ.B[50000:60000], title = "B"),
#     plot(Oₘ.t[50000:60000], Oₘ.D[50000:60000], title = "D"),
# )

X = Matrix(Oₘ[1:1000:(end-1),1:3])
dXdt = (Matrix(Oₘ[2:1000:end,1:3]) - Matrix(Oₘ[1:1000:(end-1),1:3])) / dt
# plot(sum.(abs2, eachrow(diff(dXdt, dims = 1))))

Θ = libararize(X, 4)

function STLSQ(Θ, dXdt, λ = 10^(-6))
    Ξ = Θ \ dXdt
    dim = size(Ξ, 2)
    __🚫 = 0
    
    while true
        print(".")
        🚫 = abs.(Ξ) .< λ
        Ξ[🚫] .= 0
        for j in 1:dim
            i_ = .!🚫[:, j]
            Ξ[i_, j] = Θ[:,i_] \ dXdt[:,j]
        end
        if __🚫 == 🚫 println("Stopped!"); break end # Early stopping
        __🚫 = deepcopy(🚫)
    end
    println("MSE: ", sum(abs2, Θ * Ξ - dXdt) / length(dXdt))

    return Ξ
end
Ξ = STLSQ(Θ, dXdt, 0.1)
# Ξ = STLSQ(Θ[1:1000,:], dXdt[1:1000,:], 50)

function polynomial(v, K)
    V = []

    for k ∈ 0:K
        for case = collect(multiexponents(length(v),k))
            push!(V, v .^ case)
        end
    end
    return V
end

function dynamics(Ξ)
    function f(v)
        V = prod.(polynomial(v, 4))
        return vec(V' * Ξ)
    end
    return f
end
f = dynamics(Ξ)

f(X[1,:])
dXdt[1,:]
f(X[1000,:])
dXdt[1000,:]
f(X[end,:])
dXdt[end,:]


function RK4(ODE::Function,v::Array{Float64,1},h=10^(-2))
    V1 = ODE(v)
    V2 = ODE(v .+ (h/2)*V1)
    V3 = ODE(v .+ (h/2)*V2)
    V4 = ODE(v .+ h*V3)
    return @. v + (h/6)*(V1 + 2*V2 + 2*V3 + V4)
end

function Euler(f, v, dt)
    return v + dt * f(v)
end


v0 = X[1,:]
v_ = [v0]
for t in 0:dt:2
    push!(v_, RK4(f, v_[end], dt))
end
v_
plot(reduce(hcat,v_)[2,:])
plot!(Oₘ.B[1:length(0:dt:2)])

f(v_[end])

f(X[1,:] + randn(3))
f.(Ref(X[1,:]) .+ 0.1*[randn(3) for _ in 1:100])

plot(reduce(hcat,f.(v_))[2,:])
plot!(dt\diff(Oₘ.B)[1:length(0:dt:2)])

L = CSV.read("test/lorenz.csv", DataFrame) |> Matrix
Ξ2 = STLSQ(libararize(L[1:(end-1), :],2), diff(L, dims = 1)/0.001, 0.1)
function g(v)
    V = prod.(polynomial(v, 2))
    return vec(V' * Ξ2)
end

v0 = L[1,:]
v_ = [v0]
for t in 1:100000
    push!(v_, Euler(g, v_[end], dt))
end
plot(eachrow(reduce(hcat,v_))...)


using StatsBase
using SparseArrays

function percentiler(v::AbstractVector; resolution = 10)
    V = []
    p_ = percentile.(Ref(v), range(0,100,resolution+1))
    # p_[end] += eps(p_[end])
    for k in 1:resolution
        push!(V, p_[k] .≤  v .< p_[k+1])
    end
    return sparse(reduce(hcat, V)), p_
end
foo, _ = percentiler(X[:,1])

function percentiler(X::AbstractMatrix; resolution = 10)
    ϕ_ = []
    p_ = []
    for v in eachcol(X)
        ϕ, p = percentiler(v, resolution = resolution)
        push!(ϕ_, ϕ)
        push!(p_, p)
    end
    # percentiler.(eachcol(X), resolution = resolution)
    return reduce(hcat, ϕ_), reduce(hcat, p_)
end
Φ, P = percentiler(X, resolution = 10)

using LinearAlgebra

STLSQ(Θ, dXdt, 0.1)
STLSQ(Φ, dXdt, 0.1) # |> heatmap
STLSQ([Θ Φ], dXdt, 0.1)

scatter(Oₘ.P[50000:100:60000], Oₘ.B[50000:100:60000])
scatter(Oₘ.B[50000:100:60000], Oₘ.D[50000:100:60000])
scatter(Oₘ.D[50000:100:60000], Oₘ.P[50000:100:60000])

function col_prod(A::AbstractMatrix, B::AbstractMatrix)
    C_ = []
    for a in eachcol(A)
        push!(C_, a .* B)
    end
    C = reduce(hcat, C_)
    if (typeof(A) <: AbstractSparseMatrix) ||
       (typeof(B) <: AbstractSparseMatrix)
        C = sparse(C)
    end
    return C
end
ΘΦ = col_prod(Θ, Φ)

# Ξ3 = STLSQ(ΘΦ, dXdt, 0.1)

Φ, P = percentiler(X, resolution = 10)
ΘΦ = col_prod(Θ, Φ)
@time Ξ3 = STLSQ(ΘΦ, dXdt, 0.1)
# MSE: 0.14891364899249673

function mynamics(Ξ3)
end
# baz = mynamics(Ξ3)

function baz(v)
    A_ = []
    for k in axes(P, 2)
        Aidx = findfirst(v[k] .< P[:,k])
        if (Aidx |> isnothing) || (Aidx ≤ 1) || (Aidx > size(P, 1))
            push!(A_, zeros(size(P, 1)))
        else
            push!(A_, axes(P, 1) .== (Aidx-1))
        end
    end
    pop!.(A_)
    return (col_prod(prod.(polynomial(v, 4))', vcat(A_...)') * Ξ3) |> vec
end

v = X[100,:]
dXdt[100,:]
baz(v)

plot(sum.(abs2, eachrow(dXdt) - baz.(eachrow(X))), yscale = :log10)
sum(sum.(abs2, eachrow(dXdt) - baz.(eachrow(X)))) / length(dXdt)

w = v + dt*baz(v)
u = w + dt*baz(w)
z = u + dt*baz(u)
y = z + dt*baz(z)

u - w # D, 47.25 is the problem
baz(u) - baz(w)

plot(Oₘ.P[(100*1000-100):(100*1000+100)])
plot(Oₘ.B[(100*1000-100):(100*1000+100)])
plot(Oₘ.D[(100*1000-100):(100*1000+100)])

v0 = X[100,:]
v_ = [v0]
for t in 0:dt:2
    push!(v_, Euler(baz, v_[end], dt))
end
v_

plot(reduce(hcat,v_)[2,:])
plot!(Oₘ.B[1:length(0:dt:2)])


using SparseArrays

# col_prod(col_prod(Φ[:, 1:10], Φ[:, 10 .+ (1:10)]), Φ[:, 20 .+ (1:10)]) |> sizeof
sΦ = col_prod(col_prod(Φ[:, 1:10], Φ[:, 10 .+ (1:10)]), Φ[:, 20 .+ (1:10)])
ΘsΦ = col_prod(Θ, sΦ)

Ξ4 = STLSQ(ΘsΦ, dXdt, 0.1) |> sparse

function bax(v)
    A_ = []
    for k in axes(P, 2)
        Aidx = findfirst(v[k] .< P[:,k])
        if (Aidx |> isnothing) || (Aidx ≤ 1) || (Aidx > size(P, 1))
            push!(A_, zeros(size(P, 1)))
        else
            push!(A_, axes(P, 1) .== (Aidx-1))
        end
    end
    pop!.(A_)
    polys = prod.(polynomial(v, 4))'
    location = reduce(col_prod, adjoint.(A_))
    return (col_prod(polys, location) * Ξ4) |> vec
end

# using ProgressBars
# for t0 = ProgressBar(1:1440)

t0 = 100
v = X[t0,:]
dXdt[t0,:]
bax(v)

v0 = X[t0,:]
v_ = [v0]
for t in 0:dt:2
    push!(v_, Euler(bax, v_[end], dt))
end
plot(title = "B")
plot!(Oₘ.B[(t0-1)*1000 .+ (1:length(0:dt:2))], label = "Data", color = :black, lw = 2)
plot!(stack(v_)[2,:], label = "Recovered", style = :dash, lw = 2)
# string("G:/DDM/atopy_tests/", t0, "png") |> png
# end