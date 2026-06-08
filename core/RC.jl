
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


"""'''''''''''''''''''''''''''''''''''''''''''''

        next generation reservoir computing

'''''''''''''''''''''''''''''''''''''''''''''"""

struct NGRC
    Wout::Matrix{Float64}
    λ::Float64
end
function selfprod(asdf)
    asdf_ = []
    for i in axes(asdf, 2)
        for j in axes(asdf, 2)
            if i < j continue end
            push!(asdf_, asdf[:, i] .* asdf[:, j])
        end
    end
    return stack(asdf_)
end
function ngrc(df, syms; λ = 0)
    Ssyms, Rsyms = syms
    S = dw(Matrix(df[:, Ssyms]))
    _R = Matrix(df[:, Rsyms])
    R = [ones(nrow(df)-1) [dw(_R) up(_R)] selfprod([dw(_R) up(_R)])]
    Wout = (R'R + λ*I) \ R'S
    return NGRC(Wout, λ)
end
function (s::NGRC)(X)
    X = Matrix(X)
    _X = [ones(size(X, 1)-1) [dw(X) up(X)] selfprod([dw(X) up(X)])]
    return _X * s.Wout
end


# using DataFrames, CSV, Combinatorics, LinearAlgebra

# # data = CSV.read("G:/seasurface/nino3.4/data_GLBv0.08.csv", DataFrame)

# M = reshape('a':'y', 5, 5)
# [up(M) dw(M)]


# [[up(M) dw(M)] selfprod([up(M) dw(M)])]

# M .* ('a':'e')
# ('a':'e') .* M

# vrbl = (["u", "v"], ["w", "x", "y", "z"])

# syms = vrbl
# Ssyms, Rsyms = syms
# λ = 1
# df = DataFrame([rand(100) for _ in 1:6], ["u", "v", "w", "x", "y", "z"])

# S = dw(Matrix(df[:, Ssyms]))
# _R = Matrix(df[:, Rsyms])
# R = [ones(nrow(df)-1) [dw(_R) up(_R)] selfprod([dw(_R) up(_R)])]

# Wout = (R'R + λ*I) \ R'S

# f = ngrc(df, vrbl)
# [Matrix(df[Not(1), Ssyms]) f(df[:, Rsyms])]

# lstm = CSV.read("C:/Users/rmsms/OneDrive/바탕 화면/lstm_output1_30_2023.01.01.csv", DataFrame)
# i_ = [1, 126, 251]
# j_ = 1:125:626
# j_ = [1, 126, 251, 376, 501, 626]

