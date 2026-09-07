load_packages([:StatsBase, :LinearAlgebra, :Random, :DataFrames])

"""

    rsq(y, ŷ)

Calculates the coefficient of determination (R-squared) between predicted values `ŷ` and actual values `y`.

"""
rsq(y, ŷ) = 1 - (sum(abs2, (y .- ŷ)) / sum(abs2, (y .- mean(y))))


ape(x, y) = abs.(x .- y) ./ abs.(x)
accuracy(x, y) = count(x .== y) / length(x)
ranking(x) = sortperm(sortperm(x, rev = true))
mae(x, y; kargs...) = mean(abs, x - y; kargs...)
mse(x, y; kargs...) = mean(abs2, x - y; kargs...)
rmse(x, y; kargs...) = sqrt.(mean(abs2, x - y; kargs...))
mape(x, y; kargs...) = mean(ape(x, y); kargs...)

# mse(f, data) = sum(abs2, stack(residual(f, data))) / prod(size(data[:, f.lname]))
# aic(f, data) = prod(size(data[:, f.lname]))*log(mse(f, data)) + 2nrow(f.recipe)

cov(x, y) = mean((x .- mean(x)) .* (y .- mean(y)))
cor(x, y) = cov(x, y) / (std(x) * std(y))

import Base.rand
rand(df::AbstractDataFrame; n::Integer = 1) = df[rand(1:nrow(df), n), :]

"""
    add_fold!(data::AbstractDataFrame; k = 5, seed = -1)

Adds a `fold` column to the DataFrame `data`, assigning each row to one of `k` folds.
The rows are shuffled and assigned to folds in a round-robin manner. If `seed` is non-negative,
it sets the random seed for reproducibility.

# Examples
```jldoctest
julia> X = DataFrame(a = 1:10, b = 11:20);
julia> add_fold!(X, k = 3, seed = 42)
julia> X
10×3 DataFrame
 Row │ a      b      fold
     │ Int64  Int64  Int64       
─────┼──────────────────────────────
   1 │     1     11      2
   2 │     2     12      1
   3 │     3     13      3
   4 │     4     14      2
   5 │     5     15      1
   6 │     6     16      3
   7 │     7     17      2
   8 │     8     18      1
   9 │     9     19      3
  10 │    10     20      2
```    
"""
function add_fold!(data::AbstractDataFrame; k = 5, seed = -1)
    seed ≥ 0 && Random.seed!(seed)
    data[!, :fold] .= 1 .+ mod.(shuffle(1:nrow(data)), k)
    return data
end

"""
    half(x)

    Splits the input vector `x` into two halves and returns them as a tuple.
"""
half(x) = (x[1:div(length(x), 2)], x[div(length(x), 2)+1:end])

"""
    up(array)

    Returns the array excluding the last row.
"""
up(matrix) = matrix[1:(end-1),:]

"""
    dw(array)

    Returns the array excluding the first row.
"""
dw(matrix) = matrix[2:end,:]

"""
    arglmax(x)

    Returns the indices of local maxima in the input vector `x`.
"""
function arglmax(x)
    if allequal(Float32.(x)) return [1] end
    if isempty(x) return Int64[] end
    bits = circshift(x, 1) .< x .> circshift(x, -1)
    bits[1] = all(x[1] .≥ x)
    bits[end] = all(x[end] .≥ x)
    return findall(bits)
end

"""
    arglmin(x)

    Returns the indices of local minima in the input vector `x`.
"""
function arglmin(x)
    if allequal(Float32.(x)) return [1] end
    if isempty(x) return Int64[] end
    bits = circshift(x, 1) .> x .< circshift(x, -1)
    bits[1] = all(x[1] .≤ x)
    bits[end] = all(x[end] .≤ x)
    return findall(bits)
end

"""
    dict2bifurcation(bfcn)

    Converts a dictionary of bifurcation data into two vectors: one for the keys (parameters) and one for the values (bifurcation points). The keys are repeated according to the length of their corresponding values.
"""
dict2bifurcation(bfcn) = [[fill(k, length(v)) for (k,v) in bfcn]...;], [values(bfcn)...;]

callbfcn() = Dict{Float64, Vector{Float64}}()
function callbfcn(filename)
    if isfile(filename)
        _bfcn = JLD2.load(filename)["bfcn"]
    else
        @info "Initializing bifurcation"
        _bfcn = callbfcn()
        JLD2.save(filename, bfcn = _bfcn)
    end
    return _bfcn
end

function getproxy(h, key)
    try
        return h[key]
    catch e
        key_h = collect(keys(h))
        value_h = collect(values(h))
        proxykey = argmin(abs.(key_h .- key))
        return value_h[proxykey]
    end
end

function cleanse(float::Type, bifurcation)
    bifurcation = deepcopy(deepcopy(bifurcation))
    for (k, v) in bifurcation
        bifurcation[k] = unique(float.(v))
    end
    return bifurcation
end
cleanse(bifurcation) = cleanse(Float16, bifurcation)