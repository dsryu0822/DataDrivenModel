Rsq(x, y) = 1 .- sum(abs2, (x .- y)) / sum(abs2, (x .- mean(x)))
APE(x, y) = abs.(x .- y) ./ abs.(x)
MAPE(x, y) = mean(APE(x, y))
MSE(x, y) = mean(abs2, (x .- y))
ranking(x) = sortperm(sortperm(x, rev = true))

cov(x, y) = mean((x .- mean(x)) .* (y .- mean(y)))
cor(x, y) = cov(x, y) / (std(x) * std(y))

import Base.rand
rand(df::AbstractDataFrame; n = 1) = df[rand(1:nrow(df), n), :]

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