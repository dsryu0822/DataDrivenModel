using Combinatorics, LinearAlgebra, SparseArrays, DataFrames

struct STLSQresult
    matrix::AbstractMatrix
    MSE::Float64
end
function Base.show(io::IO, ed::STLSQresult)
    show(io, "text/plain", sparse(ed.matrix))
    # println()
    print(io, "\nMSE = $(ed.MSE)")    
end

function STLSQ(ΘX, Ẋ; λ = 10^(-6), verbose = false)::STLSQresult
    if ΘX isa AbstractSparseMatrix
        ΘX = Matrix(ΘX)
    end
    Ξ = ΘX \ Ẋ
    dim = size(Ξ, 2)
    __🚫 = 0
    
    while true
        verbose && print(".")
        🚫 = abs.(Ξ) .< λ
        Ξ[🚫] .= 0
        for j in 1:dim
            i_ = .!🚫[:, j]
            Ξ[i_, j] = ΘX[:,i_] \ Ẋ[:,j]
        end
        if __🚫 == 🚫 verbose && println("Stopped!"); break end # Earl_X stopping
        __🚫 = deepcopy(🚫)
    end
    MSE = sum(abs2, Ẋ - ΘX * Ξ) / length(Ẋ)
    verbose && println("MSE = $MSE")

    return STLSQresult(Ξ, MSE)
end
function STLSQ(df::AbstractDataFrame, Ysyms::AbstractVector{T}, Xsyms::AbstractVector{T};
    K = 1, f_ = nothing,
    λ = 10^(-6), verbose = false) where T <: Union{Integer, Symbol}
    X = Θ(df[:, Xsyms], K = K, f_ = f_)
    Y = Matrix(df[:, Ysyms])
    return STLSQ(X, Y, λ = λ, verbose = verbose)
end
# X = rand(0:9, 5,3)
# v = X[1,:]
# STLSQ(X,v,0.1)

function Θ(X::AbstractMatrix; K = 1, f_ = nothing)
    dim = size(X, 2)
    ansatz = []

    for k in 0:K
        for case = collect(multiexponents(dim, k))
            push!(ansatz, prod(X .^ case', dims = 2))
        end
    end
    ΘX = hcat(ansatz...)
    if f_ |> !isnothing
        for f in f_
            ΘX = [ΘX f.(X)]
        end
    end
    return ΘX
end
Θ(X::AbstractVector; K = 1, f_ = nothing) = Θ(reshape(X, 1, :), K = K, f_ = f_)
Θ(X::AbstractDataFrame; K = 1, f_ = nothing) = Θ(Matrix(X), K = K, f_ = f_)
Θ(X::DataFrameRow; K = 1, f_ = nothing) = Θ(collect(X), K = K, f_ = f_)