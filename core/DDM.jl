struct STLSQresult
    recipe::AbstractDataFrame
    recipeF::AbstractDataFrame # fast version of `recipe`
    matrix::AbstractMatrix
    matrixF::AbstractMatrix # fast version of `matrix`
    mse::Float64
    aic::Float64
    r2::Float64
    lname::AbstractVector
    rname::AbstractVector
end
function Base.show(io::IO, s::STLSQresult)
    # show(io, "f_with_nz$(nrow(s.recipe))")
    print(io, "f_with_nz$(nrow(s.recipe))")
    # show(io, "text/plain", s.matrix)
    # print(io, "\nnumber of terms: ", nrow(s.recipe), "\t mse = ", s.mse)
end
function (s::STLSQresult)(x) # fast but unstable
    return vec(Θ(x, s.recipeF) * s.matrixF)
end
function (s::STLSQresult)(data::AbstractDataFrame)
    fitted = s(Matrix(data[:, s.rname]))
    if typeof(fitted) <: AbstractVector
        fitted = reshape(fitted, :, 1)
    end
    return DataFrame(reshape(fitted, :, length(s.lname)), s.lname)
end


function STLSQ(ΘX, Ẋ; λ = 0, verbose = false)
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
function SINDy(df::AbstractDataFrame, sysms::Tuple, recipe::AbstractDataFrame; λ = 0)
    Ysyms, Xsyms = sysms
    X = Θ(df[:, Xsyms], recipe)
    Y = Matrix(df[:, Ysyms])
    Ξ = STLSQ(X, Y, λ = λ)
    bit_sparse = all.(map(x -> iszero.(x), eachrow(Ξ)))
    # sparse_rows = findall(bit_sparse)
    recipeF = recipe[.!bit_sparse, :]
    _Ξ = Ξ[.!bit_sparse, :]
    mse = sum(abs2, Y - X * Ξ) / length(Y) # compare to original data
    aic = length(Y) * log(mse) + 2nrow(recipe)
    r2 = Rsq(vec(X * Ξ), vec(Y))
    return STLSQresult(recipe, recipeF, Ξ, _Ξ, mse, aic, r2, Ysyms, Xsyms)
end
function SINDy()
    return STLSQresult(DataFrame(), DataFrame(), zeros(1,1), zeros(1,1), 0.0, 0.0, 0.0, [], [])
end

# TODO: Add ISTA (Iterative Shrinkage-Thresholding Algorithm) for L1 regularization

# function soft_thresholding(x, λ)
#     return sign(x) * max(abs(x) - λ, 0.0)
# end
# function ISTA(A, b, λ, α=1.0, max_iter=1000)
#     x = zeros(size(A, 2))
#     for i in 1:max_iter
#         grad = A' * (A * x - b)
#         x = soft_thresholding.(x - α * grad, α * λ)
#     end
#     return x
# end


"""
    residual(f, df)

Calculate the residual of the SINDy model `f` with respect to the data frame `df`.
"""
function residual(f, df)
    # return f.(eachrow(df[:, f.rname])) - collect.(eachrow(df[:, f.lname]))
    return Matrix(f(df)) - Matrix(df[:, f.lname])
end

function pretty_term(xyz)
    xyz = replace(xyz, "sin1" => "sin", "cos1" => "cos", "pi1" => "π", "pi" => "π")
    dict_sup = Dict("0" => "⁰", "1" => "¹", "2" => "²", "3" => "³", "4" => "⁴", "5" => "⁵", "6" => "⁶", "7" => "⁷", "8" => "⁸", "9" => "⁹")
    xyz_ = split(xyz, '⋅')
    xyz1 = sort(unique(xyz_))
    deg_ = [count(xyz_ .== x) for x in xyz1]
    str_deg_ = [replace(string(d), dict_sup...) for d in deg_]
    str_deg_ = replace(str_deg_, "¹" => "")
    return join(xyz1 .* str_deg_, '⋅')
end
constant(x) = 1
for m in 1:(2^10) eval(Meta.parse("_$(m)(x) = $m * x")) end

"""
    cook(xnames; poly = 0:1, trig = 0:0, f_ = [], format = sinpi)

Generate a recipe for SINDy based on the specified polynomial and trigonometric terms.
- `xnames`: Names of the variables.
- `poly`: Degrees of polynomial terms to include.
- `trig`: Frequencies of trigonometric terms to include.
- `f_`: Custom functions to include.
- `format`: Format for trigonometric functions (e.g., `sin`, `cos`, `sinpi`, `cospi`).

"""
function cook(xnames; poly = 0:1, trig = 0:0, f_ = [], format = sinpi)
    xnames = string.(xnames)
    recipe = DataFrame(deg = Int64[], term = String[], func = [], funh = [], vecv = [])
    for d = UnitRange(extrema(poly)...)
        if d == 1
            for i in eachindex(xnames)
                push!(recipe, [1, xnames[i], identity, vec, [i]])
            end
        else
            for i in eachindex(xnames)
                rcp = recipe[recipe.deg .== (d-1), :]
                for j in i:nrow(rcp)
                    push!(recipe, [d, "$(xnames[i])⋅$(rcp.term[j])", prod, eachrow, sort([i; rcp.vecv[j]])])
                end
            end
        end
    end
    if 0 ∈ poly
        pushfirst!(recipe, [0, "1", constant, vec, [1]])
    end
    recipe = recipe[recipe.deg .∈ Ref(poly), :]

    fsin, fcos = format ∈ [sin, cos] ? (sin, cos) : format ∈ [sind, cosd] ? (sind, cosd) : (sinpi, cospi)
    for m in trig
        iszero(m) && continue
        for i in eachindex(xnames)
            push!(recipe, [m, string(fsin, m, xnames[i]), fsin, eval(Meta.parse("_$m")), [i]])
            push!(recipe, [m, string(fcos, m, xnames[i]), fcos, eval(Meta.parse("_$m")), [i]])
        end
    end
    for f in f_
        for i in eachindex(xnames)
            push!(recipe, [0, "$(string(f))($(xnames[i]))", f, vec, [i]])
        end
    end
    recipe.term = pretty_term.(recipe.term)
    unique!(recipe, :term)
    insertcols!(recipe, 1, :index => 1:nrow(recipe))

    return recipe[:, [:index, :term, :func, :funh, :vecv]]
end

Θ(X, recipe) = hcat([dr.func.(dr.funh(X[:, dr.vecv])) for dr in eachrow(recipe)]...)
Θ(X::AbstractVector, recipe) = Θ(reshape(X, 1, :), recipe)
Θ(X::AbstractDataFrame, recipe) = Θ(Matrix(X), recipe)
Θ(X::DataFrameRow, recipe) = Θ(collect(X), recipe)

import Base: print, println
import PrettyTables: pretty_table
print(s::STLSQresult) = print(pretty_table(s))
println(s::STLSQresult) = println(pretty_table(s))
function pretty_table(s::STLSQresult)
    table = [eachindex(s.recipe.term) s.recipe.term s.matrix]
    table[table .== 0] .= ""
    return pretty_table(String, table; column_labels = ["idx"; "basis"; string.(s.lname)])
end

function jacobian(T::Type, s::STLSQresult)
    rname = eval(Meta.parse("@variables $(join(string.(s.rname), " "))"))
    fnexp = vec(sum(Θ(rname, s.recipeF)' .* s.matrixF, dims = 1))

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
    mutate(s, k)

Mutates the `k`-th element of the boolean vector `s` by flipping its value.
"""
function mutate(s, k)
    _s = deepcopy(s)
    _s[k] = .!_s[k]
    return _s
end
"""
    forwardselect(DATA, vrbl, cnfg)

Performs forward selection on the datasets using the specified variables and configuration.
Returns the configuration that minimizes the Akaike Information Criterion (AIC).
"""
function forwardselect(cnfg, data, vrbl; maxitr = 1000)
    jumpt = detect_jump(data, vrbl)
    sets = set_divider(jumpt)
    datasets = [data[set, :] for set in sets]
    sampled = [ds[centerpick(nrow(ds), 50), :] for ds in datasets]

    taboo = Dict()
    strand_ = []
    for i in 1:nrow(cnfg)
        strand = zeros(Bool, nrow(cnfg))
        strand = mutate(strand, i)

        aic_ = [Inf]
        for itr in 1:maxitr
            if itr == maxitr @warn "Forward selection reached maximum iterations without convergence." end
            
            aic = fill(Inf, nrow(cnfg))
            for k in eachindex(strand)
                mutation = mutate(strand, k)
                iszero(mutation) && continue
                if hash(mutation) ∈ keys(taboo)
                    aic[k] = taboo[hash(mutation)]
                else
                    _cnfg = cnfg[mutation, :]           
                    any([rank(Θ(smpl, _cnfg)) for smpl in sampled] .< count(mutation)) && continue

                    _aic = [SINDy(smpl, vrbl, _cnfg).aic for smpl in sampled]
                    aic[k] = sum(_aic)
                    taboo[hash(mutation)] = aic[k]
                end
            end
            strand = mutate(strand, argmin(aic))
            if minimum(aic) < minimum(aic_)
                push!(aic_, minimum(aic))
            else
                break
            end
        end
        push!(strand_, minimum(aic_) => strand)
    end
    strand = last.(strand_)[argmin(first.(strand_))]
    
    return cnfg[strand, :]
end

"""
    centerpick(n, m)

Select `m` evenly spaced indices from `1:n`, excluding the first and last indices.
"""
centerpick(n, m) = round.(Int64, range(1, n, m+2))[2:end-1]
"""
    labeling!(data, vrbl, cnfg; θ = 0)

Label the subsystems in DataFrame `data` with respect to `vrbl` and `cnfg` configuration.
`dos` is the degree of smoothness.
"""
function labeling!(data, vrbl, cnfg; θ = 0)
    jumpt = detect_jump(data, vrbl)
    sets = set_divider(jumpt)
    datasets = [data[set, :] for set in sets]
    sampled = [ds[centerpick(nrow(ds), 50), :] for ds in datasets]
    # cnfg = forwardselect(sampled, vrbl, cnfg)

    f_ = [SINDy(smpl, vrbl, cnfg) for smpl in sampled]
    L = zeros(length(f_), length(f_))
    for i in eachindex(f_)
        for j in eachindex(f_)
            if i > j
                L[i, j] = sum(abs2, f_[i].matrix - f_[j].matrix)
                L[j, i] = L[i, j]
            end
        end
    end
    # scatter(filter(!iszero, vec(L)), yscale = :log10)
    subsys = connected_components(SimpleGraph(L .< θ))
    f_ = [SINDy([sampled[ss]...;], vrbl, cnfg) for ss in subsys]; # print.(f_)

    label = zeros(Int64, nrow(data))
    for i in eachindex(subsys)
        for j in subsys[i]
            label[sets[j]] .= i
        end
    end

    bit_blank = iszero.(label)
    label[bit_blank] .= argmin.(eachrow(sum.(abs2, stack([residual(f, data[bit_blank, :]) for f in f_]))))
    data.label = label

    return data
end


function dryad(data, cols; method = :tree) # fairy of tree and forest
    # _data = deepcopy(data[(data.label .!= circshift(data.label, -1)) .|| (data.label .!= circshift(data.label, 1)), :])
    # add_fold!(_data)
    apply_ = method == :forest ? apply_forest : apply_tree
    build_ = method == :forest ? build_forest : build_tree

    labels = data.label
    features = Matrix(data[:, cols])
    acc_ = []
    for seed in 1:10
        Dtree = build_(data.label, features, rng = seed); # print_tree(Dtree, feature_names = ["V", "I", "Vr"])
        push!(acc_, count(labels .== apply_(Dtree, features)) / length(labels))
        if maximum(acc_) ≈ 1 break end #; else print("█") end
    end
    Dtree = build_(data.label, features, rng = argmax(acc_))
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



function export_tex(s::STLSQresult; digits = 4)
    coeflines = vec.(eachcol(s.matrixF))
    bases = s.recipeF.term
    texstring = ""
    for (coefline, LHS) in zip(coeflines, s.lname)
        texstring *= "$LHS =&"
        for (coef, basis) in zip(coefline, bases)
            if iszero(coef)
                continue
            else
                texstring *= " + $(round(coef, digits=digits))\\,$basis"
            end
        end
        texstring *= " \\\\\n"
        texstring = replace(texstring, "& +" => "&")
        texstring = replace(texstring, "+ -" => "- ")
    end
    texstring = replace(texstring, Dict(["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"] .=> ["^{$k}" for k in 0:9])...)
    return texstring
end

function export_py(s::STLSQresult; digits = 99)
    coeflines = vec.(eachcol(s.sparse_matrix))
    bases = Θ(string.(s.rname), N = s.N, M = s.M, f_ = s.f_, C = s.C)
    texstring = "def f($(join(string.(s.rname), ", "))):\n\t return ["
    for (coefline, LHS) = zip(coeflines, s.lname)
        texstring *= ""
        for (coef, basis) = zip(coefline, bases)
            if coef |> iszero
                continue
            else
                texstring *= " + $(round.(coef; digits))*$(replace(basis, " " => "*"))"
            end
        end
        texstring *= ", "
        texstring = replace(texstring, "+ -" => "- ")
    end
    texstring = replace(texstring, "[ + " => "[")
    texstring = replace(texstring, "[ - " => "[")
    texstring = replace(texstring * "]", ", ]" => "]")
    texstring = replace(texstring, Dict(["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"] .=> ["**$k" for k in 0:9])...)
    return texstring
end