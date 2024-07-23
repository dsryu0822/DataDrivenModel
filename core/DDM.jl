using Combinatorics, LinearAlgebra, SparseArrays, DataFrames, PrettyTables, Symbolics
# @info "Combinatorics, LinearAlgebra, SparseArrays, DataFrames, PrettyTables, Symbolics loaded"

struct STLSQresult
    N::Int64
    M::Int64
    f_::Array{Function}
    matrix::AbstractMatrix
    MSE::Float64
    lname::AbstractVector
    rname::AbstractVector
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
    Ξ = sparse(Ξ ./ L₂) # L₂ is row-wise producted to denormalize coefficient matrix
    return Ξ
end
function SINDy(X::AbstractMatrix, Ẋ::AbstractMatrix;
    N = 1, M = 0, f_ = Function[],
    λ = 10^(-6), verbose = false)

    ΘX = Θ(X; N = N, M = M, f_ = f_)
    Ξ = STLSQ(ΘX, Ẋ, λ = λ, verbose = verbose)
    MSE = sum(abs2, Ẋ - _ΘX * Ξ) / length(Ẋ) # compare to original data
    lname = "dx" .* string.(axes(Ẋ, 2))
    rname =  "x" .* string.(axes(X, 2))
    return STLSQresult(N, M, f_, Ξ, MSE, lname, rname)
end
function SINDy(df::AbstractDataFrame, Ysyms::AbstractVector{T}, Xsyms::AbstractVector{T};
    N = 1, M = 0, f_ = Function[],
    λ = 10^(-6), verbose = false) where T <: Union{Integer, Symbol}

    X = Θ(df[:, Xsyms], N = N, M = M, f_ = f_)
    Y = Matrix(df[:, Ysyms])
    Ξ = STLSQ(X, Y, λ = λ, verbose = verbose)
    MSE = sum(abs2, Y - X * Ξ) / length(Y) # compare to original data
    return STLSQresult(N, M, f_, Ξ, MSE, Ysyms, Xsyms)
end


function Θ(X::AbstractMatrix; N = 1, M = 0, f_ = Function[], λ = 0)
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

const supdigit = Dict(0:9 .=> ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"])
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
function print(s::STLSQresult)
    table = [1:size(s.matrix, 1) Θ(string.(s.rname), N = s.N, M = s.M, f_ = s.f_) s.matrix]
    table[table .== 0] .= ""
    return pretty_table(table; header = ["idx"; "basis"; string.(s.lname)])
end

function jacobian(s::STLSQresult)
    # lname = eval(Meta.parse("@variables $(join(string.(s.lname), " "))"))
    rname = eval(Meta.parse("@variables $(join(string.(s.rname), " "))"))
    fnexp = vec(sum(Θ(rname, N = s.N, M = s.M, f_ = s.f_)' .* s.matrix, dims = 1))
    return Symbolics.jacobian(fnexp, rname)
end


# function FDM1(M::AbstractMatrix, dt = 0.1)
#     d = size(M, 2)
#     names = [fill("x", d) .* string.(1:d); fill("dx", d) .* string.(1:d)]
#     return DataFrame([diff(M, dims = 1)/dt M[2:end, :]], names)
# end

function add_subsystem!(data, vrbl, cnfg; θ1 = 1e-1, θ2 = 1e-24, θ3 = 1e-10, min_rank = 0)
    normeddf = sum.(abs2, eachrow(diff(Matrix(data[:, first(vrbl)]), dims = 1))) # scatter(normeddf[1:100:end], yscale = :log10)
    jumpt = [1; findall(normeddf .> θ1)]
    _sets = filter(!isempty, collect.(UnitRange.(jumpt .+ 1, circshift(jumpt .- 1, -1)))); pop!(_sets); _sets = UnitRange.(first.(_sets), last.(_sets))
    sets = filter(!isempty, sort(union(_sets, UnitRange.(last.(_sets)[1:(end-1)] .+ 2, first.(_sets)[2:end] .- 2))))

    subsystem = zeros(Int64, nrow(data));
    for id_subsys = 1:8 # id_subsys = 1; id_subsys = 2; id_subsys = 3
        rank_ = [rank(Θ(Matrix(data[a_set,last(vrbl)]); cnfg...)) for a_set in sets]
        if maximum(rank_) < min_rank
            # candy = SINDy(data[iszero.(subsystem),:], vrbl...; cnfg...); print(candy, last(vrbl))
            for (A, B) = combinations(sets, 2)
                candy = SINDy([data[A, :]; data[B, :]], vrbl...; cnfg...)
                # A = first(sets); B = last(sets);
                if candy.MSE < θ2 break end
            end
            # @warn "No pair found!"
            # if candy.MSE ≥ θ2
            for (A, B, C) = combinations(sets, 3)
                candy = SINDy([data[A, :]; data[B, :]; data[C, :]], vrbl...; cnfg...)
                # A = first(sets); B = last(sets);
                if candy.MSE < θ2 break end
            end
            # end
        else
            meaningful = sets[argmax(rank_)]
            candy = SINDy(data[meaningful,:], vrbl...; cnfg...)
        end
        # print(candy)

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