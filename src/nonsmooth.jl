using Combinatorics
using StatsBase
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

X = rand(0:9, 5,3)
v = X[1,:]

function poly_basis(X::AbstractMatrix, K = 1)
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
poly_basis(X, 4)

function poly_basis(v::AbstractVector, K = 1)
    ansatz = []

    for k ∈ 0:K
        for case = collect(multiexponents(length(v),k))
            push!(ansatz, v .^ case)
        end
    end
    return vcat(ansatz...)
end
poly_basis(v, 4)

# function poly_dynamics(Ξ)
#     function f(v)
#         V = prod.(polynomial(v, 4))
#         return vec(V' * Ξ)
#     end
#     return f
# end


function RK4(f::Function,v::Array{Float64,1}, h=10^(-2))
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4)
end

function Euler(f::Function,v::Array{Float64,1}, h=10^(-2))
    return v + h*f(v)
end

function percentiler(v::AbstractVector; resolution = 10)
    V = []
    p_ = percentile.(Ref(v), range(0,100,resolution+1))
    p_[end] += eps(p_[end])
    for k in 1:resolution
        push!(V, p_[k] .≤  v .< p_[k+1])
    end
    return sparse(reduce(hcat, V)), p_
end
function percentiler(X::AbstractMatrix; resolution = 10)
    ϕ_ = []
    p_ = []
    for v in eachcol(X)
        ϕ, p = percentiler(v, resolution = resolution)
        push!(ϕ_, ϕ)
        push!(p_, p)
    end
    return ϕ_, p_
    # return reduce(hcat, ϕ_), reduce(hcat, p_)
end

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
    return (col_prod(polys, location) * Ξ) |> vec
end
