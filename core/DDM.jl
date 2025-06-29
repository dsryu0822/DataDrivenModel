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
(s::STLSQresult)(data::DataFrameRow) = s(collect(data[s.rname]))
(s::STLSQresult)(data::AbstractDataFrame) = vcat(s.(eachrow(data[:, s.rname]))...)


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

"""
    residual(f, df)

Calculate the residual of the SINDy model `f` with respect to the data frame `df`.
"""
function residual(f, df)
    return sum.(abs2, f.(eachrow(df[:, f.rname])) - collect.(eachrow(df[:, f.lname])))
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
    nr, nc = size(X)
    padding = zeros(nr)
    nz_ = Int64[]
    i = 0
    incest = UnitRange.(1, cumsum([length(multiexponents(nc+1, N)); length(f_)*nc; 2M*nc]))
    incest = [Int64.(incest[1]), setdiff(incest[2], incest[1]), setdiff(incest[3], incest[2])]

    ansatz = []
    for k in 0:N
        for case = collect(multiexponents(nc, k))
            i += 1
            if i ∈ sparse_rows θx = padding else
                push!(nz_, i)
                θx = prod(X .^ case', dims = 2)
            end
            push!(ansatz, θx)
        end
    end
    ΘX = hcat(ansatz...)
    for x ∈ eachcol(X)
        for f in f_
            i += 1
            if i ∈ sparse_rows θx = padding else
                push!(nz_, i)
                θx = f.(x)
            end
            ΘX = [ΘX θx]
        end
        for m in 1:M
            i += 1
            if i ∈ sparse_rows θx = padding else
                push!(nz_, i)
                θx = cospi.(m*x)
            end
            ΘX = [ΘX θx]
        end
        for m in 1:M
            i += 1
            if i ∈ sparse_rows θx = padding else
                push!(nz_, i)
                θx = sinpi.(m*x)
            end
            ΘX = [ΘX θx]
        end
    end

    for c in 2:C
        for j_ in combinations(2:size(ΘX, 2), c)
            if !any(Ref(j_) .⊆ incest)
                i += 1
                if i ∈ sparse_rows θx = padding else
                    push!(nz_, i)
                    θx = .*([ΘX[:, j] for j in j_]...)
                end
                ΘX = [ΘX θx]
            end
        end
    end

    return ΘX[:, nz_]
end
   Θ(X::AbstractVector; kargs...) = Θ(reshape(X, 1, :); kargs...)
Θ(X::AbstractDataFrame; kargs...) = Θ(Matrix(X); kargs...)
     Θ(X::DataFrameRow; kargs...) = Θ(collect(X); kargs...)
function Θ(X::Vector{String}; N = 1, M = 0, f_ = Function[], C = 1, λ = 0)
    # λ is just for dummy argument for add_subsystem! function
    dim = length(X)
    ΘX = []
    incest = UnitRange.(1, cumsum([length(multiexponents(dim+1, N)); length(f_)*dim; 2M*dim]))
    incest = [Int64.(incest[1]), setdiff(incest[2], incest[1]), setdiff(incest[3], incest[2])]

    for k in 0:N
        for case = collect(multiexponents(dim, k))
            # push!(ΘX, reduce(*, ((X .* num2sup.(case))[.!iszero.(case)])))
            push!(ΘX, join(((X .* num2sup.(case))[.!iszero.(case)]), " "))
        end
    end
    for x in X
        for f in f_
            push!(ΘX, "$(string(f))($x)")
        end
        for m in 1:M
            _m = ifelse(m |> isone, "", string(m))
            push!(ΘX, ("cos$(_m)π$x"))
        end
        for m in 1:M
            _m = ifelse(m |> isone, "", string(m))
            push!(ΘX, ("sin$(_m)π$x"))
        end
    end

    dim = length(ΘX)
    for c in 2:C
        for j_ in combinations(2:dim, c)
            if !any(Ref(j_) .⊆ incest)
                push!(ΘX, *([ΘX[j] for j in j_]...))
            end
        end
    end

    replace!(ΘX, "" => "1")
    return ΘX
end
import Base: print, println
import PrettyTables: pretty_table
print(s::STLSQresult) = print(pretty_table(s))
println(s::STLSQresult) = println(pretty_table(s))
function pretty_table(s::STLSQresult)
    table = [1:size(s.sparse_matrix, 1) Θ(string.(s.rname), N = s.N, M = s.M, f_ = s.f_, C = s.C) s.sparse_matrix]
    table[table .== 0] .= ""
    return pretty_table(String, table; header = ["idx"; "basis"; string.(s.lname)])
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
    # return sets[sortperm(length.(sets), rev=true)]
    return sets
end
function detect_jump(data, vrbl; dos = 0)
# if dos == 0
    normeddf = norm.(eachrow(diff(Matrix(data[:, first(vrbl)]), dims = 1)))
# elseif dos == 1
#     normeddf = norm.(eachrow(diff(diff(Matrix(data[:, first(vrbl)]), dims = 1), dims = 1))) # scatter(normeddf[1:100:end], yscale = :log10)
# end
    len_normeddf = length(normeddf)
    _jumpt = [-1]; jumpt = deepcopy(_jumpt);
    # if jumpt is initialized with [0], then it will be a problem when idx = 1
    for _ in 1:nrow(data)
        jumpt = deepcopy(_jumpt)
        idx = argmax(normeddf)
        normeddf[idx] = -Inf

        idx3 = [max(1, idx-1), idx, min(idx+1, len_normeddf)] # min(idx+1, len_normeddf) is to prevent BoundsError
        if all(abs.(_jumpt .- idx) .> 1) && isempty(idx3 ∩ argmax(normeddf))
            push!(_jumpt, idx3...)
        else
            break
        end
    end
    jumpt = unique([1; (sort(_jumpt[2:end])); nrow(data)])
    
    # normeddf = norm.(eachrow(diff(Matrix(data[:, first(vrbl)]), dims = 1)))
    # plot(yscale = :log10, msw = 0, legend = :none);
    # # plot(yscale = :log10, msw = 0, xlims = [0, 10], legend = :none);
    # scatter!(normeddf, shape = :pixel);
    # scatter!(jumpt[1:end-1], normeddf[jumpt[1:end-1]], shape = :x);
    # png("normeddf")
    return jumpt
end

"""
    labeling!(data, vrbl, cnfg; dos = 0)

Label the subsystems in DataFrame `data` with respect to `vrbl` and `cnfg` configuration.
`dos` is the degree of smoothness.
"""
function labeling!(data, vrbl, cnfg; θ = 1e-24, dos = 0)
    jumpt = detect_jump(data, vrbl; dos)
    sets = set_divider(jumpt)
    datasets = [data[set, :] for set in sets]

    n = length(datasets)
    min_m = SINDy(datasets[1][1:10, :], vrbl...; cnfg...).sparse_matrix.m
    sampled = [ds[centerpick(nrow(ds), min_m), :] for ds in datasets]
    f__ = fill(SINDy(datasets[1][1:10, :], vrbl...; cnfg...), n, n)
    mse__ = fill(Inf, n, n)
    dist__ = zeros(n, n)
    for (i,j) = doublerange(length(datasets))
        if i ≤ j continue end
        f__[i, j] = SINDy([sampled[i]; sampled[j]], vrbl...; cnfg...)
        f__[j, i] = f__[i, j]
        mse__[i, j] = f__[i, j].MSE
        mse__[j, i] = mse__[i, j]
    end
    subsys = connected_components(SimpleGraph(mse__ .< θ))
    f_ = [SINDy([sampled[ss]...;], vrbl...; cnfg...) for ss in subsys]; # print.(f_)

    label = zeros(Int64, nrow(data))
    for i in eachindex(subsys)
        for j in subsys[i]
            label[sets[j]] .= i
        end
    end
    bit_blank = iszero.(label)
    label[bit_blank] .= argmin.(eachrow(stack([residual(f, data[bit_blank, :]) for f in f_])))
    data.label = label

    return data
end

function dryad(data, vrbl) # fairy of tree and forest
    labels = data.label
    features = Matrix(data[:, vrbl])
    acc_ = []
    for seed in 1:10
        Dtree = build_tree(data.label, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
        push!(acc_, count(labels .== apply_tree(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break end #; else print("█") end
    end
    Dtree = build_tree(data.label, features, rng = argmax(acc_))
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
