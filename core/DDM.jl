using Combinatorics, LinearAlgebra, SparseArrays, DataFrames, PrettyTables, Symbolics
# @info "Combinatorics, LinearAlgebra, SparseArrays, DataFrames, PrettyTables, Symbolics loaded"

struct STLSQresult
    N::Int64
    M::Int64
    f_::Array{Function}
    C::Int64
    matrix::AbstractMatrix
    MSE::Float64
    lname::AbstractVector
    rname::AbstractVector
end
function Base.show(io::IO, s::STLSQresult)
    show(io, "text/plain", sparse(s.matrix))
    # println()
    print(io, "\npoly N = $(s.N), fourier M = $(s.M), f_ = $(s.f_), C = $(s.C) with MSE = $(s.MSE)")
    print(io, "\nplease use print function to show all result")
end
function (s::STLSQresult)(x)
    return vec(Θ(x; N = s.N, M = s.M, f_ = s.f_, C = s.C) * s.matrix)
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
function SINDy(X::AbstractMatrix, Ẋ::AbstractMatrix;
    λ = 1e-6, verbose = false, N = 1, M = 0, f_ = [], C = 0)

    ΘX = Θ(X; N = N, M = M, f_ = f_, C = C)
    Ξ = STLSQ(ΘX, Ẋ, λ = λ, verbose = verbose)
    MSE = sum(abs2, Ẋ - ΘX * Ξ) / length(Ẋ) # compare to original data
    lname = "dx" .* string.(axes(Ẋ, 2))
    rname =  "x" .* string.(axes(X, 2))
    return STLSQresult(N, M, f_, C, Ξ, MSE, lname, rname)
end
function SINDy(df::AbstractDataFrame, Ysyms::AbstractVector{T}, Xsyms::AbstractVector{T};
    λ = 1e-6, verbose = false, N = 1, M = 0, f_ = [], C = 0) where T <: Union{Integer, Symbol}

    X = Θ(df[:, Xsyms]; N = N, M = M, f_ = f_, C = C)
    Y = Matrix(df[:, Ysyms])
    Ξ = STLSQ(X, Y, λ = λ, verbose = verbose)
    MSE = sum(abs2, Y - X * Ξ) / length(Y) # compare to original data
    return STLSQresult(N, M, f_, C, Ξ, MSE, Ysyms, Xsyms)
end


const superdigit = Dict(0:9 .=> ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"])
function num2sup(num)
    if (num == 0) || (num == 1)
        return ""
    else
        return reduce(*, (getindex.(Ref(superdigit), reverse(digits(num, base = 10)))))
    end
end
function Θ(X::AbstractMatrix; N = 1, M = 0, f_ = Function[], C = 1, λ = 0)
    # λ is just for dummy argument for add_subsystem! function
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

    dim = size(ΘX, 2)
    for c in 2:C
        for (j1, j2) in combinations(2:dim, c)
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
    table = [1:size(s.matrix, 1) Θ(string.(s.rname), N = s.N, M = s.M, f_ = s.f_, C = s.C) s.matrix]
    table[table .== 0] .= ""
    return pretty_table(table; header = ["idx"; "basis"; string.(s.lname)])
end

function jacobian(s::STLSQresult)
    # lname = eval(Meta.parse("@variables $(join(string.(s.lname), " "))"))
    rname = eval(Meta.parse("@variables $(join(string.(s.rname), " "))"))
    fnexp = vec(sum(Θ(rname, N = s.N, M = s.M, f_ = s.f_, C = s.C)' .* s.matrix, dims = 1))
    return Symbolics.jacobian(fnexp, rname)
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
function add_subsystem!(data, vrbl, cnfg; θ1 = 1e-1, θ2 = 1e-24, θ3 = 1e-10, min_rank = 0, dos = 0)
    if dos == 0
        normeddf = sum.(abs2, eachrow(diff(Matrix(data[:, first(vrbl)]), dims = 1)))
        # jumpt = [1; findall(normeddf .> θ1)]
        # _sets = filter(!isempty, collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1)))); pop!(_sets); _sets = UnitRange.(first.(_sets), last.(_sets))
        # sets = filter(!isempty, sort(union(_sets, UnitRange.(last.(_sets)[1:(end-1)] .+ 2, first.(_sets)[2:end] .- 2))))
    elseif dos == 1
        normeddf = sum.(abs2, eachrow(diff(diff(Matrix(data[:, first(vrbl)]), dims = 1), dims = 1))) # scatter(normeddf[1:100:end], yscale = :log10)
    end
    jumpt = [0]
    while true
        idx = argmax(normeddf)
        if all(abs.(jumpt .- idx) .> 1)
            push!(jumpt, idx, idx+1, idx-1)
            normeddf[[idx, idx+1, idx-1]] .= -Inf
        else
            break
        end
    end
    jumpt = [1; sort(jumpt[2:end]); nrow(data)]
    sets = set_divider(jumpt)

    subsystem = zeros(Int64, nrow(data));
    for id_subsys = 1:8 # id_subsys = 1; id_subsys = 2; id_subsys = 3
        rank_ = [rank(Θ(Matrix(data[a_set,last(vrbl)]); cnfg...)) for a_set in sets]
        if maximum(rank_) < min_rank
            for (A, B) = combinations(sets, 2)
                # A = sets[1]; B = sets[2];
                candy = SINDy([data[A, :]; data[B, :]], vrbl...; cnfg...)
                if candy.MSE < θ2 break end
            end
            for (A, B, C) = combinations(sets, 3)
                candy = SINDy([data[A, :]; data[B, :]; data[C, :]], vrbl...; cnfg...)
                if candy.MSE < θ2 break end
            end
        else
            meaningful = sets[argmax(rank_)]
            candy = SINDy(data[meaningful,:], vrbl...; cnfg...)
            # print(candy)
        end

        idx_blank = findall(iszero.(subsystem))
        residual = norm.(eachrow(Matrix(data[idx_blank, first(vrbl)])) .- candy.(eachrow(Matrix(data[idx_blank, last(vrbl)]))))
        # scatter(residual)
        # scatter(residual[1:100:end], yscale = :log10)
        # kmeaned = kmeans(reshape(log10.(residual), 1, :), 2)
        # θ3 = exp10(mean(kmeaned.centers))
        idx_blank = idx_blank[residual .< θ3]
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