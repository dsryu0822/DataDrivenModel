load_packages([:DataFrames, :LinearAlgebra])

"""'''''''''''''''''''''''''''''''''''''''''''''

                reservoir computing

'''''''''''''''''''''''''''''''''''''''''''''"""

struct Reservoir
    M::Int64
    P::Int64
    N::Int64
    α::Float64
    β::Float64
    ρ::Float64
    A::AbstractMatrix
    ξ::Float64
    Win::AbstractMatrix
    Wout::AbstractMatrix
    c::AbstractVector
    r::AbstractVector
end
function Base.show(io::IO, s::Reservoir)
    Base.print(io, "f: R^$(s.M) → R^$(s.P)")
end
function Base.print(s::Reservoir)
    Base.print("""
    f: R^$(s.M) → R^$(s.P)
      reservoir size N = $(s.N)
       leakage reate α = $(s.α)
    sparsity of Wout β = $(s.β)
     spectral radius ρ = $(s.ρ)
    """)
    return nothing
end
function (s::Reservoir)(U)    
    encode(r, u) = (1-s.α)*r + s.α*tanh.(s.A * r + s.Win * u .+ s.ξ)

    r_ = [s.r]
    [push!(r_, encode(r_[end], u)) for u in eachcol(U)]
    popfirst!(r_); R = stack(r_)
    S_pred = s.Wout * R .+ s.c
    return S_pred    
end
function (s::Reservoir)(U::DataFrame)
    return s(Matrix(U)')
end

function reservoir_computing(U::AbstractMatrix, S::AbstractMatrix;
    seed = 0, warmup = 10,
    N = 500, D = 2, α = 0.5, β = 1e-2, ρ = 1.0, σ = 1.0, ξ = 1.0)

    if seed ≥ 0 Random.seed!(seed) end
    M = size(U, 1)
    P = size(S, 1)
    RN = erdos_renyi(N, (D*N ÷ 2))
    A = adjacency_matrix(RN) .* 2(rand(N, N) .- .5)
    A = ρ*A / maximum(abs.(eigen(Matrix(A)).values))
    Win = σ*(sparse(stack([shuffle([1; zeros(M-1)]) for _ in 1:N])') .* 2(rand(N, M) .- .5))
    encode(r, u) = (1-α)*r + α*tanh.(A * r + Win * u .+ ξ)
    r_ = [zeros(N)]
    [push!(r_, encode(r_[end], u)) for u in eachcol(U)]
    popfirst!(r_); R = stack(r_)
    R = R[:, (warmup+1):end]
    S = S[:, (warmup+1):end]
    Wout = S*R'inv(R*R' + (β*LinearAlgebra.I))
    c = -vec(Wout * mean(R, dims = 2) - mean(S, dims = 2))
    return Reservoir(M, P, N, α, β, ρ, A, ξ, Win, Wout, c, r_[end])
end
function reservoir_computing(df, Sidcs, Uidcs; kargs...)
    S = Matrix(df[:, Sidcs])'
    U = Matrix(df[:, Uidcs])'
    return reservoir_computing(U, S; kargs...)
end


"""'''''''''''''''''''''''''''''''''''''''''''''

        next generation reservoir computing

'''''''''''''''''''''''''''''''''''''''''''''"""
struct NGRC
    Wout::Matrix{Float64}
    λ::Float64
end
function Base.show(io::IO, s::NGRC)
    Base.print(io, "f: R^$(size(s.Wout, 1)) → R^$(size(s.Wout, 2))")
end
function selfprod(M)
    m = size(M, 2)
    P = similar(M, eltype(M), size(M, 1), m * (m + 1) ÷ 2)
    k = 1
    @inbounds for i in 1:m
        for j in 1:m
            if i < j continue end
            P[:, k] = M[:, i] .* M[:, j]
            k += 1
        end
    end
    return P
end
function tikhonov_λ(X, y)
    m, n = size(X)
    if m < n
        F = svd(X')
        U, s, V = F.V, F.S, F.U
    else
        U, s, V = svd(X)
    end
    Uᵀy = U'y
    return λ -> V * (Diagonal(s ./ (s.^2 .+ λ)) * Uᵀy)
end
function ngrc_sweep(df, Sidcs, Uidcs)
    X = Matrix(df)
    S = dw(X[:, Sidcs]) # target
    U = X[:, Uidcs]     # observer
    R = selfprod([ones(size(U, 1)-1) dw(U) up(U)])
    # Wout = (R'R + λ*I) \ R'S
    # return NGRC(Wout, λ)
    return tikhonov_λ(R, S)
end
function ngrc(df, Sindcs, Uincs; λ = 0)
    X = Matrix(df)
    S = dw(X[:, Sindcs]) # target
    U = X[:, Uincs]     # observer
    R = selfprod([ones(size(U, 1)-1) dw(U) up(U)])
    Wout = (R'R + λ*I) \ R'S
    return NGRC(Wout, λ)
end
function (s::NGRC)(X)
    U = Matrix(X)
    R = selfprod([ones(size(U, 1)-1) dw(U) up(U)])
    return R * s.Wout
end