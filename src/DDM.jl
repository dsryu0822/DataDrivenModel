using Combinatorics
using LinearAlgebra
using SparseArrays

function STLSQ(Θ, dXdt, λ = 10^(-6))
    if Θ isa AbstractSparseMatrix
        Θ = Matrix(Θ)
    end
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

function col_prod(A::AbstractMatrix, B::AbstractMatrix)
    C_ = []
    for a in eachcol(A)
        push!(C_, a .* B)
    end
    C = reduce(hcat, C_)
    if (A isa AbstractSparseMatrix) ||
       (B isa AbstractSparseMatrix)
        C = sparse(C)
    end
    return C
end
col_prod(A::AbstractVector, B::AbstractVector) = col_prod(A', B')