using Combinatorics, LinearAlgebra, SparseArrays
using Printf

struct STLSQresult
    matrix::AbstractMatrix
    MSE::Float64
end
function Base.show(io::IO, ed::STLSQresult)
    show(io, "text/plain", sparse(ed.matrix))
    println()
    print(io, "MSE = $(ed.MSE)")    
end

function STLSQ(Θ, dXdt; λ = 10^(-6), verbose = false)::STLSQresult
    if Θ isa AbstractSparseMatrix
        Θ = Matrix(Θ)
    end
    Ξ = Θ \ dXdt
    dim = size(Ξ, 2)
    __🚫 = 0
    
    while true
        verbose && print(".")
        🚫 = abs.(Ξ) .< λ
        Ξ[🚫] .= 0
        for j in 1:dim
            i_ = .!🚫[:, j]
            Ξ[i_, j] = Θ[:,i_] \ dXdt[:,j]
        end
        if __🚫 == 🚫 verbose && println("Stopped!"); break end # Earl_X stopping
        __🚫 = deepcopy(🚫)
    end
    MSE = sum(abs2, Θ * Ξ - dXdt) / length(dXdt)
    verbose && println("MSE: $MSE")

    return STLSQresult(Ξ, MSE)
end
function STLSQ(df::AbstractDataFrame, Ysyms::AbstractVector{T}, Xsyms::AbstractVector{T},
    K = 1, f_ = nothing;
    λ = 10^(-6), verbose = false) where T <: Union{Integer, Symbol}
    if f_ |> isnothing
        X = poly_basis(Matrix(df[:, Xsyms]), K)
    else
        X = col_func(Matrix(df[:, Xsyms]), f_)
    end
    Y = Matrix(df[:, Ysyms])
    return STLSQ(X, Y, λ = λ, verbose = verbose)
end
# X = rand(0:9, 5,3)
# v = X[1,:]
# STLSQ(X,v,0.1)

function poly_basis(X::AbstractMatrix, K = 1; forcing = false)
    mtrx_X = Matrix(X)
    if size(X, 2) > size(X, 1) && !forcing
        @warn "X is not a tall matrix!"
        return nothing
    end
    dim = size(X, 2)
    ansatz = []

    for k in 0:K
        for case = collect(multiexponents(dim,k))
            push!(ansatz, prod(mtrx_X .^ case', dims = 2))
        end
    end
    return hcat(ansatz...)
end
# poly_basis(X, 4)

function poly_basis(v::AbstractVector, K = 1; forcing = false)
    if length(v) > 10 && !forcing
        @warn "X seems to be no vector!"
        return nothing
    end
    ansatz = []

    for k ∈ 0:K
        for case = collect(multiexponents(length(v),k))
            push!(ansatz, v .^ case)
        end
    end
    return prod.(ansatz)
end
# poly_basis(v, 4)

# function col_prod(A::AbstractMatrix, B::AbstractMatrix)
#     C_ = []
#     for a in eachcol(A)
#         push!(C_, a .* B)
#     end
#     C = reduce(hcat, C_)
#     if (A isa AbstractSparseMatrix) ||
#        (B isa AbstractSparseMatrix)
#         C = sparse(C)
#     end
#     return C
# end
# col_prod(A::AbstractVector, B::AbstractVector) = col_prod(A', B')

function col_func(X::AbstractMatrix, f_)
    _X = deepcopy(X)
    for f in f_
        _X = [_X f.(X)]
    end
    return _X
end

function col_func(X::AbstractVector, f_)
    _X = deepcopy(X)
    for f in f_
        _X = [_X; f.(X)]
    end
    return _X
end