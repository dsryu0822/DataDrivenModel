using Combinatorics, LinearAlgebra, SparseArrays, DataFrames, PrettyTables
@info "Packages Combinatorics, LinearAlgebra, SparseArrays
      DataFrames, PrettyTables loaded"

struct STLSQresult
    N::Int64
    M::Int64
    f_::Array{Function}
    matrix::AbstractMatrix
    MSE::Float64
end
function Base.show(io::IO, s::STLSQresult)
    show(io, "text/plain", sparse(s.matrix))
    # println()
    print(io, "\npoly N = $(s.N), fourier M = $(s.M), f_ = $(s.f_) with MSE = $(s.MSE)")
    print(io, "\nplease use print function to show all result")
end
function (s::STLSQresult)(x)
    return vec(Θ(x; N = s.N, M = s.M, f_ = s.f_) * s.matrix)
end

function STLSQ(ΘX, Ẋ; λ = 10^(-6), verbose = false)
    _ΘX = deepcopy(ΘX)
    L₂ = norm.(eachcol(ΘX))
    ΘX = ΘX ./ L₂'
    # L₂ is for column-wise normalization to ensure restricted isometry property
    # Due to this L₂, λ thresholding would be doesn't work as expected

    Ξ = ΘX \ Ẋ; dim = size(Ξ, 2)
    _🚫 = 0
    while true
        verbose && print(".")
        🚫 = abs.(Ξ) .< (λ * L₂)
        Ξ[🚫] .= 0
        for j in 1:dim
            i_ = .!🚫[:, j]
            Ξ[i_, j] = ΘX[:,i_] \ Ẋ[:,j]
        end
        if _🚫 == 🚫 verbose && println("Stopped!"); break end # Earl_X stopping
        _🚫 = deepcopy(🚫)
    end
    Ξ =  sparse(Ξ ./ L₂) # L₂ is row-wise producted to denormalize coefficient matrix
    MSE = sum(abs2, Ẋ - _ΘX * Ξ) / length(Ẋ) # compare to original data
    verbose && println("MSE = $MSE")

    return Ξ, MSE
end
function SINDy(X::AbstractMatrix, Ẋ::AbstractMatrix;
    N = 1, M = 0, f_ = Function[],
    λ = 10^(-6), verbose = false)
    ΘX = Θ(X; N = N, M = M, f_ = f_)
    Ξ, MSE = STLSQ(ΘX, Ẋ, λ = λ, verbose = verbose)
    return STLSQresult(N, M, f_, Ξ, MSE)
end
function SINDy(df::AbstractDataFrame, Ysyms::AbstractVector{T}, Xsyms::AbstractVector{T};
    N = 1, M = 0, f_ = Function[],
    λ = 10^(-6), verbose = false) where T <: Union{Integer, Symbol}
    X = Θ(df[:, Xsyms], N = N, M = M, f_ = f_)
    Y = Matrix(df[:, Ysyms])
    Ξ, MSE = STLSQ(X, Y, λ = λ, verbose = verbose)
    return STLSQresult(N, M, f_, Ξ, MSE)
end
# X = rand(0:9, 5,3)
# v = X[1,:]
# STLSQ(X,v,0.1)


function Θ(X::AbstractMatrix; N = 1, M = 0, f_ = Function[])
    dim = size(X, 2)
    ansatz = []

    for k in 0:N
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
   Θ(X::AbstractVector; N = 1, M = 0, f_ = Function[]) = 
    Θ(reshape(X, 1, :), N = N, M = M, f_ = f_)
Θ(X::AbstractDataFrame; N = 1, M = 0, f_ = Function[]) = 
           Θ(Matrix(X), N = N, M = M, f_ = f_)
     Θ(X::DataFrameRow; N = 1, M = 0, f_ = Function[]) = 
          Θ(collect(X), N = N, M = M, f_ = f_)

supdigit = Dict(0:9 .=> ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"])
function num2sup(num)
    if (num == 0) || (num == 1)
        return ""
    else
        return reduce(*, (getindex.(Ref(supdigit), reverse(digits(num, base = 10)))))
    end
end

function Θ(X::Vector{String}; N = 1, M = 0, f_ = Function[])
    dim = length(X)
    ΘX = []

    for k in 0:N
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
function print(s::STLSQresult, vals::AbstractVector{String})
    @assert !isempty(vals) "empty vals"
    table = [1:size(s.matrix, 1) Θ(vals, N = s.N, M = s.M, f_ = s.f_) s.matrix]
    table[table .== 0] .= ""
    return pretty_table(table; header = ["idx"; "basis"; "d" .* vals[1:size(s.matrix, 2)]])
end
print(S::STLSQresult, vals::AbstractVector{Symbol}) = print(S, string.(vals))

function FDM1(M::AbstractMatrix, dt = 0.1)
    d = size(M, 2)
    names = [fill("x", d) .* string.(1:d); fill("dx", d) .* string.(1:d)]
    return DataFrame([diff(M, dims = 1)/dt M[2:end, :]], names)
end

# ets = SINDy(rand(10, 3), rand(10, 3), N = 2, M = 3, f_ = [tan])
# print(ets, ["u", "v", "x"])
# ets([1,2,3])