using Combinatorics, LinearAlgebra, SparseArrays, DataFrames, PrettyTables
@info "Packages Combinatorics, LinearAlgebra, SparseArrays
      DataFrames, PrettyTables loaded"

struct STLSQresult
    K::Int64
    M::Int64
    f_::Array{Function}
    matrix::AbstractMatrix
    MSE::Float64
end
function Base.show(io::IO, s::STLSQresult)
    show(io, "text/plain", sparse(s.matrix))
    # println()
    print(io, "\npoly K = $(s.K), fourier M = $(s.M), f_ = $(s.f_) with MSE = $(s.MSE)")
    print(io, "\nplease use print function to show all result")
end
function (s::STLSQresult)(x)
    return vec(Θ(x; K = s.K, M = s.M, f_ = s.f_) * s.matrix)
end

function STLSQ(ΘX, Ẋ; λ = 10^(-6), verbose = false)
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
    Ξ =  sparse(Ξ)
    MSE = sum(abs2, Ẋ - ΘX * Ξ) / length(Ẋ)
    verbose && println("MSE = $MSE")

    return Ξ, MSE
end
function SINDy(X::AbstractMatrix, Ẋ::AbstractMatrix;
    K = 1, M = 0, f_ = Function[],
    λ = 10^(-6), verbose = false)
    ΘX = Θ(X; K = K, M = M, f_ = f_)
    Ξ, MSE = STLSQ(ΘX, Ẋ, λ = λ, verbose = verbose)
    return STLSQresult(K, M, f_, Ξ, MSE)
end
function SINDy(df::AbstractDataFrame, Ysyms::AbstractVector{T}, Xsyms::AbstractVector{T};
    K = 1, M = 0, f_ = Function[],
    λ = 10^(-6), verbose = false) where T <: Union{Integer, Symbol}
    X = Θ(df[:, Xsyms], K = K, M = M, f_ = f_)
    Y = Matrix(df[:, Ysyms])
    Ξ, MSE = STLSQ(X, Y, λ = λ, verbose = verbose)
    return STLSQresult(K, M, f_, Ξ, MSE)
end
# X = rand(0:9, 5,3)
# v = X[1,:]
# STLSQ(X,v,0.1)


function Θ(X::AbstractMatrix; K = 1, M = 0, f_ = Function[])
    dim = size(X, 2)
    ansatz = []

    for k in 0:K
        for case = collect(multiexponents(dim, k))
            push!(ansatz, prod(X .^ case', dims = 2))
        end
    end
    ΘX = hcat(ansatz...)
    for f in f_
        ΘX = [ΘX f.(X)]
    end
    for m in 1:M
        ΘX = [ΘX cospi.(m*X) sinpi.(m*X)]
    end

    return ΘX
end
   Θ(X::AbstractVector; K = 1, M = 0, f_ = Function[]) = 
    Θ(reshape(X, 1, :), K = K, M = M, f_ = f_)
Θ(X::AbstractDataFrame; K = 1, M = 0, f_ = Function[]) = 
           Θ(Matrix(X), K = K, M = M, f_ = f_)
     Θ(X::DataFrameRow; K = 1, M = 0, f_ = Function[]) = 
          Θ(collect(X), K = K, M = M, f_ = f_)

supdigit = Dict(0:9 .=> ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"])
function num2sup(num)
    if (num == 0) || (num == 1)
        return ""
    else
        return reduce(*, (getindex.(Ref(supdigit), reverse(digits(num, base = 10)))))
    end
end

function Θ(X::Vector{String}; K = 1, M = 0, f_ = Function[])
    dim = length(X)
    ΘX = []

    for k in 0:K
        for case = collect(multiexponents(dim, k))
            push!(ΘX, reduce(*, ((X .* num2sup.(case))[.!iszero.(case)])))
        end
    end
    for f in f_
        push!(ΘX, (string(f) .* "(" .* X .* ")")...)
    end
    for m in 1:M
        _m = ifelse(m |> isone, "", string(m))
        push!(ΘX, ("cos$(_m)π" .* X)..., ("sin$(_m)π" .* X)...)
    end

    # ΘX = lpad.(ΘX, maximum(length.(ΘX)))
    replace!(ΘX, "" => "1")
    return ΘX
end

import Base: print
function print(s::STLSQresult, vals)
    @assert !isempty(vals) "empty vals"
    table = [1:size(s.matrix, 1) Θ(vals, K = s.K, M = s.M, f_ = s.f_) s.matrix]
    table[table .== 0] .= ""
    return pretty_table(table; header = ["idx"; "basis"; "d" .* vals])
end

function FDM1(M::AbstractMatrix, dt = 0.1)
    d = size(M, 2)
    names = [fill("x", d) .* string.(1:d); fill("dx", d) .* string.(1:d)]
    return DataFrame([diff(M, dims = 1)/dt M[2:end, :]], names)
end

# ets = SINDy(rand(10, 3), rand(10, 3), K = 2, M = 3, f_ = [tan])
# print(ets, ["u", "v", "x"])
# ets([1,2,3])