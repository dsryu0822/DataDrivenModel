include("../src/factorio.jl")
include("../src/DDM.jl")
include("../src/visual.jl")

function arg1(df)
    x = sum.(abs2, eachrow(df))
    y = abs.((x ./ circshift(x, 1)) .- 1)
    y[1] = minimum(y[2:end])
    return y
end

function argjump(df, thld = -8)
    y = arg1(df)
    z = findall(log10.(y) .> thld)
    return z, plot(log10.(y))
end

function maxdiff(z)
    tk = argmax(diff(z))
    return z[tk], z[tk+1]
end
# dfdistance(A, B, cols) = (Matrix(A[:, cols]) - Matrix(B[:, cols])) |> eachrow .|> norm |> maximum
dfdistance(A, B, cols) = (Matrix(A[:, cols]) - Matrix(B[:, cols])) |> eachrow .|> norm
timemean(y) = cumsum(y) ./ (1:length(y))

using LinearAlgebra, DecisionTree

d_range = 0.1:0.0001:0.3
plan = DataFrame(idx=eachindex(d_range), d=d_range)

cd("//155.230.155.221/ty/DS");
pwd()

dr = eachrow(plan)[1]
@time data = factory_soft(dr.idx, dr.d)
@time data = factory_soft(dr.idx, dr.d, (45, 60))
test = 500000

s_ = unique(round.(Int64, 10 .^ (0:0.1:4)))
err1_ = []
err2_ = []
부검 = []
# sparsity = 1
@time for sparsity = ProgressBar(s_)
_data = data[1:sparsity:test, :]
tend = nrow(_data)
# tend = nrow(_data)
# x = maxdiff(argjump(_data, -3))
# sampled = rand(UnitRange(x...), 20)
# sampled = UnitRange(x...)
sampled = rand(1:tend, 20)
STLSQed = STLSQ(_data[sampled, :], [:du, :dv], [:t, :u, :v], M = 2, verbose=true, λ = 0.001)
while STLSQed.MSE > 0.0001
    sampled = rand(1:tend, 20)
    STLSQed = STLSQ(_data[sampled, :], [:du, :dv], [:t, :u, :v], M = 2, verbose=true, λ = 0.001)
end
α = STLSQed.matrix

residual = Matrix(_data[:, [:du, :dv]]) - (Θ(_data[:, [:t, :u, :v]], M = 2) * α) |>
           eachrow .|> norm .|> log10
_data[:, :subsystem] = 1 .+ (residual .> -5)
scatter(residual[1:10:end])

STLSQ_ = [STLSQ(gdf, [:du, :dv], [:t, :u, :v], M = 2, λ = 0.001)
          for gdf in groupby(_data, :subsystem)]

# push!(err1_, getproperty.(STLSQ_, :MSE) |> sum)
# push!(부검, STLSQ_)

function factory_STLSQ(STLSQed)
    function f(s, x)
        return [1; vec(Θ(x, M = 2) * STLSQed[s].matrix)]
    end
    return f
end
g = factory_STLSQ(STLSQ_)

dtree = DecisionTreeClassifier(max_depth=3)
fit!(dtree, Matrix(_data[:, [:t, :u, :v]]), _data.subsystem)
print_tree(dtree)

x_ = [collect(data[test+1, [:t, :u, :v]])]
x = x_[end]

try
    # for t in ProgressBar(x[1]:10^(-5):55)
    for t in x[1]:10^(-5):60
        s = predict(dtree, x_[end])
        x, dx = RK4(g, s, x_[end], 10^(-5))
        push!(x_, x)
    end
    pred = DataFrame(stack(x_)', [:t, :u, :v])
    push!(err2_, dfdistance(data[(test + 1):end, :], pred, [:u, :v]))
catch e
    push!(err2_, NaN)
    # @info "DomainError for sparsity = $sparsity"
end
end

plot(log10.(last.(timemean.(err2_))))
# plot(findfirst.(map(x -> x .> 0.6, err2_)) ./ 500000,
# xlabel = "Sparsity " * L"S", title = "First passage: error > 0.5")

# plot(timemean(dfdistance(data[(test + 1):end, :], pred, [:u, :v])))

# plot(cumsum(dfdistance(data[(test + 1):end, :], pred, [:u, :v]))[1:100:end])
# plot(data[(test+1):10:end, :u], data[(test+1):10:end, :v])
# plot!(pred[1:10:end, :u], pred[1:10:end, :v], color = :red)
# plot(Matrix(data[(test+1):end, [:u, :v]]) - Matrix(pred[:, [:u, :v]]) |> eachrow .|> norm)
# plot()
# maximum(err2_[end])
# findfirst(err2_[end] .> .63)

# plot(err2_[end][1:10:end])
# plot(err2_[end][1:10:end], xlims = (11000, 12000))

# plot(data[(test + 1):10:end, :].v, xlims = (11000, 12000))
# plot!(pred[1:10:end, :].v, color = :red, xlims = (11000, 12000))


# plot(data[(test + 1):10:end, :].v)
# plot!(pred[1:10:end, :].v, color = :red)

pargs = (xscale = :log10, xlims = (1, 10^4), legend = :none, yscale = :log10)
p1 = plot(s_, err1_; pargs...)
title!(p1, "Residual")

p2 = scatter(s_, last.(timemean.(err2_)); pargs...)
xlabel!(p2, "Sparsity " * L"S")
ylabel!(p2, "Error")
vline!(p2, [500], color = :red)
# title!(p2, "Maximal error " * L"|| x - \hat{x} || _\infty")

plot(p1, p2, layout = (2,1))
