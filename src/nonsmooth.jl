using Combinatorics
using SparseArrays
using LinearAlgebra

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


function argjump(type::Type, Y::AbstractVector)
    Δ²Y = diff(diff(Y))
    Δ²Y = Δ²Y ./ maximum(Δ²Y)
    bit_jump = Δ²Y .> 0.1
    println(sum(bit_jump), " jumps detected!")
    if type == Bool
        return bit_jump
    elseif type == Int64
        return findall(bit_jump)
    else
        @error "type must be Bool or Int64"
    end
end
argjump(Y) = argjump(Int64, Y) # If type is not specified, return indices

# using StatsBase
# function percentizer(v::AbstractVector; nclass = 10)
#     V = []
#     p_ = percentile.(Ref(v), range(0,100,nclass+1))
#     p_[end] += eps(p_[end])
#     for k in 1:nclass
#         push!(V, p_[k] .≤  v .< p_[k+1])
#     end
# return sparse(reduce(hcat, V)), p_
# end
# """
# nclass is a number or sorted float array.
# """
# function percentizer(X::AbstractMatrix; nclass = 10)
#     ϕ_ = []
#     p_ = []
#     for v in eachcol(X)
#         ϕ, p = gridize(v, nclass = nclass)
#         push!(ϕ_, ϕ)
#         push!(p_, p)
#     end
#     return ϕ_, p_
#     # return reduce(hcat, ϕ_), reduce(hcat, p_)
# end
# function bax(v)
#     A_ = []
#     for k in axes(P, 2)
#         Aidx = findfirst(v[k] .< P[:,k])
#         if (Aidx |> isnothing) || (Aidx ≤ 1) || (Aidx > size(P, 1))
#             push!(A_, zeros(size(P, 1)))
#         else
#             push!(A_, axes(P, 1) .== (Aidx-1))
#         end
#     end
#     pop!.(A_)
#     polys = prod.(polynomial(v, 4))'
#     location = reduce(col_prod, adjoint.(A_))
#     return (col_prod(polys, location) * Ξ) |> vec
# end
