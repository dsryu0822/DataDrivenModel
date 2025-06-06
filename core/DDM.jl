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
    for f in f_
        push!(ΘX, (string(f) .* "(" .* X .* ")")...)
    end
    for m in 1:M
        _m = ifelse(m |> isone, "", string(m))
        push!(ΘX, ("cos$(_m)π" .* X)..., ("sin$(_m)π" .* X)...)
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
    if dos == 0
        normeddf = norm.(eachrow(diff(Matrix(data[:, first(vrbl)]), dims = 1)))
    elseif dos == 1
        normeddf = norm.(eachrow(diff(diff(Matrix(data[:, first(vrbl)]), dims = 1), dims = 1))) # scatter(normeddf[1:100:end], yscale = :log10)
    end
    len_normeddf = length(normeddf)
    _jumpt = [-1]; jumpt = deepcopy(_jumpt);
    # if jumpt is initialized with [0], then it will be a problem when idx = 1
    while true
        idx = argmax(normeddf)
        if all(abs.(_jumpt .- idx) .> 1)
            jumpt = deepcopy(_jumpt)
            idx3 = [max(1, idx-1), idx, min(idx+1, len_normeddf)]
            # min(idx+1, len_normeddf) is to prevent BoundsError
            push!(_jumpt, idx3...)
            normeddf[idx3] .= -Inf
            # println(idx3)
        else
            break
        end
    end
    jumpt = unique([1; (sort(_jumpt[2:end])); nrow(data)])
    
    # normeddf = norm.(eachrow(diff(diff(Matrix(data[:, first(vrbl)]), dims = 1), dims = 1)))
    # plot(yscale = :log10, msw = 0, legend = :none);
    # # plot(yscale = :log10, msw = 0, xlims = [0, 10], legend = :none);
    # scatter!(normeddf, shape = :+);
    # scatter!(jumpt[1:end], normeddf[jumpt[1:end-1]], shape = :x);
    # png("normeddf")
    return jumpt
end

# """
#     add_subsystem!(data, vrbl, cnfg; θ = 1e-24, dos = 0)

# Add subsystem to DataFrame `data` with respect to `vrbl` and `cnfg` configuration.
# `θ` is the threshold for residual error and `dos` is the degree of smoothness.
# """
# function add_subsystem!(data, vrbl, cnfg; θ = 1e-24, dos = 0)
#     jumpt = detect_jump(data, vrbl; dos)

#     sets = set_divider(jumpt)
#     subsystem = zeros(Int64, nrow(data));
#     # candy = SINDy(data[last(sets), :], vrbl...; cnfg...)
#     candy = nothing
#     for id_subsys = 1:6 # id_subsys = 0; id_subsys += 1
#         flag = false
        
#         # if subsystem |> iszero
#         #     candy = SINDy(data[first(sets), :], vrbl...; cnfg...)
#         # else
#         #     candy = SINDy(data[iszero.(subsystem), :], vrbl...; cnfg...)
#         # end
#         # if candy.MSE > θ
#             for many = 3:-1:1 # many = 1; many = 2; many = 3; cane = first(combinations(sets, many))
#                 for cane = combinations(sets, many)
#                     sugar = reduce(vcat, [data[cn, :] for cn in cane])
#                     candy = SINDy(sugar, vrbl...; cnfg...)
#                     if candy.MSE < θ
#                         flag = true
#                         break
#                     end
#                 end
#                 if flag break end
#             end
#         # end

#         idx_blank = findall(iszero.(subsystem))
#         residual = sum.(abs2, eachrow(Matrix(data[idx_blank, first(vrbl)])) .- candy.(eachrow(Matrix(data[idx_blank, last(vrbl)]))))
#         # scatter(residual[1:100:end], yscale = :log10); hline!([θ], color = :red)
#         idx_blank = idx_blank[residual .< θ]
#         subsystem[idx_blank] .= id_subsys
#         sets = sets[getindex.(sets, length.(sets) .÷ 2) .∉ Ref(idx_blank)]
#         # candy
        
#         if sets |> isempty break end
#     end
#     data[!, :subsystem] = subsystem;
#     return data
# end

"""
    labeling!(data, vrbl, cnfg; dos = 0)

Label the subsystems in DataFrame `data` with respect to `vrbl` and `cnfg` configuration.
`dos` is the degree of smoothness.
"""
function labeling!(data, vrbl, cnfg; dos = 0)
    jumpt = detect_jump(data, vrbl; dos = dos)

    sets = set_divider(jumpt)
    datasets = [data[set, :] for set in sets]
ground_truth = [ifelse(rand(ds.u) > 0.05, 2, ifelse(rand(ds.u) < -0.05, 3, 1)) for ds in datasets]    

    # subsystem = zeros(Int64, nrow(data));
    label = zeros(Int64, length(sets))
    idcs = eachindex(sets)
    
    for id_subsys in eachindex(sets)
        idx_hero = argmax(iszero.(label) .* length.(sets))
        aliens = setdiff(idcs, [idx_hero] .+ [-1, 0, 1])

        report_hero = DataFrame(alien = Int64[], mse = Float64[])
        for alien = aliens
            data_alien = [datasets[[alien, idx_hero]]...;]
            f = SINDy(data_alien, vrbl...; cnfg...)
            push!(report_hero, [alien, f.MSE])
        end
rargs = (; xticks = 1:30, xlims = [0, 30], yscale = :log10, ylabel = "MSE", legend = :none, ms = 12, msw = 0, ma = .5,)
s1 = scatter(report_hero.alien, report_hero.mse, text = ground_truth[report_hero.alien]; rargs...)
if id_subsys == 1 png(s1, "hyperparameter/details/$cnfg-itr_$(id_subsys).png") end

        idx_hive = report_hero.alien[argmax(report_hero.mse)]
        report_hive = DataFrame(alien = Int64[], f = [], mse = Float64[])
        for alien = aliens
            # if idx_hive == alien continue end
            data_hive = [datasets[[alien, idx_hive]]...;]
            f = SINDy(data_hive, vrbl...; cnfg...)
            push!(report_hive, [alien, f, f.MSE])
        end
# s2 = scatter(s1, report_hive.alien, report_hive.mse, text = ground_truth[report_hive.alien]; rargs...)
        worst = SINDy([datasets[idcs]...;], vrbl...; cnfg...)
        if maximum(report_hero.mse) < worst.MSE
            idx_team = report_hero.mse .< report_hive.mse
        else
            idx_team = report_hero.mse .≥ 0
        end
# s3 = scatter(s1, report_hero.alien[idx_team], report_hero.mse[idx_team], shape = :x; rargs...)

        idx_labeled = [idx_hero; report_hero.alien[idx_team]]
        label[idx_labeled] .= id_subsys
        idcs = setdiff(idcs, idx_labeled)
        if length(idcs) == 0
            break
        elseif length(idcs) == 1
            label[only(idcs)] = id_subsys + 1
            break
        end
    end

    data[!, :label] = zeros(Int64, nrow(data))
    for (k, lbl) in enumerate(label)
        data.label[sets[k]] .= lbl
    end
    bit_zero = iszero.(data.label)
    f_ = [SINDy(data, vrbl...; cnfg...) for data in groupby(data[.!bit_zero,:], :label)]
    data.label[bit_zero] .= argmin.(eachrow(stack([sum.(abs2, f.(eachrow(data[bit_zero, last(vrbl)])) - collect.(eachrow(data[bit_zero, first(vrbl)]))) for f in f_])))

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
