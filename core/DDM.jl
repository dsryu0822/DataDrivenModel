using Combinatorics, LinearAlgebra, SparseArrays, DataFrames, PrettyTables, Symbolics
# @info "Combinatorics, LinearAlgebra, SparseArrays, DataFrames, PrettyTables, Symbolics loaded"

struct STLSQresult
    N::Int64
    M::Int64
    f_::Array{Function}
    C::Int64
    sparse_matrix::AbstractMatrix
    sparse_rows::AbstractVector
    dense_matrix::AbstractMatrix
    MSE::Float64
    lname::AbstractVector
    rname::AbstractVector
end
function Base.show(io::IO, s::STLSQresult)
    show(io, "text/plain", s.sparse_matrix)
    # println()
    print(io, "\npoly N = $(s.N), fourier M = $(s.M), f_ = $(s.f_), C = $(s.C) with MSE = $(s.MSE)")
    print(io, "\nplease use print function to show all result")
end
# function (s::STLSQresult)(x) # slow but stable and concrete
#     return vec(Θ(x; N = s.N, M = s.M, f_ = s.f_, C = s.C) * s.sparse_matrix)
# end
function (s::STLSQresult)(x) # fast but unstable
    return vec(Θ(x; N = s.N, M = s.M, f_ = s.f_, C = s.C, sparse_rows = s.sparse_rows) * s.dense_matrix)
end

# function functionalizer(s::STLSQresult) # x4 slower than direct matrix multiplication
#     rname = eval(Meta.parse("@variables $(join(string.(s.rname), " "))"))
#     fnexp = vec(sum(Θ(rname, N = s.N, M = s.M, f_ = s.f_, C = s.C)' .* s.matrix, dims = 1))
#     return v -> substitute(fnexp, Dict(rname .=> v))
# end

function STLSQ(ΘX, Ẋ; λ = 1e-6, verbose = false)
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
        if _🚫 == 🚫 verbose && println("Stopped!"); break end
        _🚫 = deepcopy(🚫)
    end
    Ξ = sparse(Ξ ./ L₂) # L₂ is row-wise producted to denormalize coefficient matrix
    return Ξ
end
function SINDy(df::AbstractDataFrame, Ysyms::AbstractVector{T}, Xsyms::AbstractVector{T};
    λ = 1e-6, verbose = false, N = 1, M = 0, f_ = [], C = 0) where T <: Union{Integer, Symbol}

    X = Θ(df[:, Xsyms]; N = N, M = M, f_ = f_, C = C)
    Y = Matrix(df[:, Ysyms])
    Ξ = STLSQ(X, Y, λ = λ, verbose = verbose)
    sparse_rows = findall(all.(map(x -> iszero.(x), eachrow(Ξ))))
    _Ξ = Ξ[.!all.(map(x -> iszero.(x), eachrow(Ξ))), :]
    MSE = sum(abs2, Y - X * Ξ) / length(Y) # compare to original data
    return STLSQresult(N, M, f_, C, Ξ, sparse_rows, _Ξ, MSE, Ysyms, Xsyms)
end


const dict_superdigit = Dict(0:9 .=> ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"])
function num2sup(num)
    if (num == 0) || (num == 1)
        return ""
    else
        return reduce(*, (getindex.(Ref(dict_superdigit), reverse(digits(num, base = 10)))))
    end
end
function Θ(X::AbstractMatrix;
    N = 1, M = 0, f_ = Function[], C = 1, λ = 0, sparse_rows = [])
    # λ is just for dummy argument for add_subsystem! function
    dim = size(X, 2)
    ansatz = []
    i = 0

    for k in 0:N
        for case = collect(multiexponents(dim, k))
            i += 1; i ∈ sparse_rows && continue
            push!(ansatz, prod(X .^ case', dims = 2))
        end
    end
    ΘX = hcat(ansatz...)
    for x ∈ eachcol(X)
        for f in f_
            i += 1; i ∈ sparse_rows && continue
            ΘX = [ΘX f.(x)]
        end
        for m in 1:M
            i += 1; i ∈ sparse_rows && continue
            ΘX = [ΘX cospi.(m*x)]
        end
        for m in 1:M
            i += 1; i ∈ sparse_rows && continue
            ΘX = [ΘX sinpi.(m*x)]
        end
    end

    dim = size(ΘX, 2)
    for c in 2:C
        for (j1, j2) in combinations(2:dim, c)
            i += 1; i ∈ sparse_rows && continue
            ΘX = [ΘX (ΘX[:, j1] .* ΘX[:, j2])]
        end
    end

    return ΘX
end
   Θ(X::AbstractVector; kargs...) = Θ(reshape(X, 1, :); kargs...)
Θ(X::AbstractDataFrame; kargs...) = Θ(Matrix(X); kargs...)
     Θ(X::DataFrameRow; kargs...) = Θ(collect(X); kargs...)
function Θ(X::Vector{String}; N = 1, M = 0, f_ = Function[], C = 1, λ = 0)
    # λ is just for dummy argument for add_subsystem! function
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

    dim = length(ΘX)
    for c in 2:C
        for (j1, j2) in combinations(2:dim, c)
            push!(ΘX, ΘX[j1] * ΘX[j2])
        end
    end

    replace!(ΘX, "" => "1")
    return ΘX
end
import Base: print
function print(s::STLSQresult)
    table = [1:size(s.sparse_matrix, 1) Θ(string.(s.rname), N = s.N, M = s.M, f_ = s.f_, C = s.C) s.sparse_matrix]
    table[table .== 0] .= ""
    return pretty_table(table; header = ["idx"; "basis"; string.(s.lname)])
end

function jacobian(T::Type, s::STLSQresult)
    rname = eval(Meta.parse("@variables $(join(string.(s.rname), " "))"))
    fnexp = vec(sum(Θ(rname, N = s.N, M = s.M, f_ = s.f_, C = s.C)' .* s.sparse_matrix, dims = 1))

    J = Symbolics.jacobian(fnexp, rname)
    if T == Matrix
        return J
    elseif T == Function
        return x -> Float64.(substitute(J, Dict(rname .=> x)))
    else
        error("Type not supported: Only `Function` or `Matrix` are supported.")
    end
end

function set_divider(arr::AbstractVector)
    s = 0
    sets = []
    flag_record = false
    for k in first(arr):last(arr)
        if flag_record
            if (k+1) ∈ arr
                push!(sets, s:k)
                flag_record = false
            end
        else
            if k ∉ arr
                s = k
                flag_record = true
            end
        end
    end
    return sets[sortperm(length.(sets), rev=true)]
end

"""
    add_subsystem!(data, vrbl, cnfg; θ = 1e-24, dos = 0)

Add subsystem to DataFrame `data` with respect to `vrbl` and `cnfg` configuration.
`θ` is the threshold for residual error and `dos` is the degree of smoothness.
"""
function add_subsystem!(data, vrbl, cnfg; θ = 1e-24, dos = 0)
    if dos == 0
        normeddf = norm.(eachrow(diff(Matrix(data[:, first(vrbl)]), dims = 1)))
        # jumpt = [1; findall(normeddf .> θ1)]
        # _sets = filter(!isempty, collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1)))); pop!(_sets); _sets = UnitRange.(first.(_sets), last.(_sets))
        # sets = filter(!isempty, sort(union(_sets, UnitRange.(last.(_sets)[1:(end-1)] .+ 2, first.(_sets)[2:end] .- 2))))
    elseif dos == 1
        normeddf = norm.(eachrow(diff(diff(Matrix(data[:, first(vrbl)]), dims = 1), dims = 1))) # scatter(normeddf[1:100:end], yscale = :log10)
    end
    _jumpt = [0]; jumpt = [0];
    while true
        idx = argmax(normeddf)
        if all(abs.(_jumpt .- idx) .> 1)
            jumpt = deepcopy(_jumpt)
            push!(_jumpt, idx, idx+1, idx-1)
            normeddf[[idx, idx+1, idx-1]] .= -Inf
        else
            break
        end
    end
    jumpt = [1; sort(jumpt[2:end]); nrow(data)]
    sets = set_divider(jumpt)

    subsystem = zeros(Int64, nrow(data));
    for id_subsys = 1:8 # id_subsys = 1; id_subsys = 2; id_subsys = 3
        flag = false
        
        candy = SINDy(data[rand(sets), :], vrbl...; cnfg...)
        for many = 1:3 # many = 1; many = 2; many = 3; cane = first(combinations(sets, many))
            for cane = combinations(sets, many)
                sugar = reduce(vcat, [data[cn, :] for cn in cane])
                candy = SINDy(sugar, vrbl...; cnfg...)
                if candy.MSE < θ
                    flag = true
                    break
                end
            end
            if flag break end
        end

        idx_blank = findall(iszero.(subsystem))
        residual = sum.(abs2, eachrow(Matrix(data[idx_blank, first(vrbl)])) .- candy.(eachrow(Matrix(data[idx_blank, last(vrbl)]))))
        # scatter(residual[1:100:end], yscale = :log10)
        idx_blank = idx_blank[residual .< θ]
        subsystem[idx_blank] .= id_subsys
        sets = sets[rand.(sets) .∉ Ref(idx_blank)]
        
        if sets |> isempty break end
        # if sets |> isempty @info "sets exhausted"; break end
    end
    data[!, :subsystem] = subsystem;
    return data
end
function dryad(data, vrbl) # fairy of tree and forest
    labels = data.subsystem
    features = Matrix(data[:, vrbl])
    acc_ = []
    for seed in 1:10
        Dtree = build_tree(data.subsystem, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
        push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break end #; else print("█") end
    end
    Dtree = build_tree(data.subsystem, features, rng = argmax(acc_))
    # println("Accuracy: $(count(labels .== apply_tree(Dtree, features)) / length(labels))")
    return Dtree
end

function gram_schmidt(J)
    N = size(J, 1)
    U, V = deepcopy(J), deepcopy(J)
    U[:,1] = V[:,1] / norm(V[:,1])
    for j in 2:N
        for jp in 1:(j-1)
            V[:,j] -= (J[:,j]'U[:,jp])*U[:,jp]
        end
        U[:,j] = V[:,j] / norm(V[:,j])
    end
    return U, V
end