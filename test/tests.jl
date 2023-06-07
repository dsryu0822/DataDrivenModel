@time using DataDrivenDiffEq
@time using ModelingToolkit
@time using OrdinaryDiffEq
@time using DataDrivenSparse
@time using LinearAlgebra
@time using Plots

# Create a test problem
function lorenz(u, p, t)
    x, y, z = u

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z
    return [ẋ, ẏ, ż]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
dt = 0.1
prob = ODEProblem(lorenz, u0, tspan)
# sol = solve(prob, Tsit5(), dense = true) # sol(sol.t,Val{4}); sol(sol.t,Val{5})
sol = solve(prob, RK4(), saveat = dt) # Tsit5() is more efficient https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Explicit-Runge-Kutta-Methods
plot(
    plot(sol, idxs = (0, 1), color = :black),
    plot(sol, idxs = (0, 2), color = :black),
    plot(sol, idxs = (0, 3), color = :black),
    plot(sol, idxs = (1, 2, 3), color = :black)
)
# How get the derivative from ODEsolution
# https://github.com/SciML/DataDrivenDiffEq.jl/blob/8acc287468ef13be134ebd2d85d273bc0f456cfb/src/problem/type.jl#L486
@assert sol(sol.t,Val{1})[1] == ((sol(dt) - sol(0)) / dt)

## Start the automatic discovery
ddprob = DataDrivenProblem(sol)

@variables x y z
u = [x; y; z]
basis = Basis(polynomial_basis(u, 2), u)
# @variables t x(t) y(t) z(t)
# u = [x; y; z]
# basis = Basis(polynomial_basis(u, 2), u, iv = t)

opt = STLSQ(exp10.(-5:0.1:-1))
# https://github.com/SciML/DataDrivenDiffEq.jl/blob/cae0c79e0d7b9eef48f9350e01f69b0798b447bb/lib/DataDrivenSparse/src/algorithms/STLSQ.jl#L97-L120
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 4))
soleq = get_basis(ddsol)
println(soleq) # soleq[1].rhs
get_parameter_map(soleq) # Equlvalent to get_parameter_values(soleq)
ddsol.residuals # Equlvalent to rss(ddsol)

@time using DataDrivenDiffEq
@time using DataDrivenSparse

@time using CSV, DataFrames, Dates
function fillmissing!(data)
    for col in eachcol(data)
        while true
            bit_missing = ismissing.(col)
            if sum(bit_missing) == 0 break end
            col[bit_missing] .= ((circshift(col, -1) + circshift(col, 1)) / 2)[bit_missing]
            bit_missing = ismissing.(col)
            col[bit_missing] .= circshift(col, 1)[bit_missing]
        end
    end
end

data = DataFrame()
DATA_ = []
for fn = readdir("data")
    itemname = Symbol(first(split(fn, " ")))
    if (itemname == :알루미늄) || (itemname == :오렌지) continue end
    push!(DATA_, CSV.read("data/" * fn, DataFrame))
    DATA_[end].날짜 = Date.(DATA_[end].날짜, dateformat"y- m- d")
    select!(DATA_[end], ["날짜", "종가"])
    rename!(DATA_[end], [:t, itemname])
    if eltype(DATA_[end][:, 2]) <: AbstractString
        DATA_[end][:, 2] = replace.(DATA_[end][:, 2], "," => "")
        DATA_[end][!, 2] = parse.(Float64, DATA_[end][:, 2])
    end
    if isempty(data)
        data = deepcopy(DATA_[end])
    else
        leftjoin!(data, DATA_[end], on = :t, makeunique = true)
    end
end
data = data[data.t .< Date(2022, 1, 1), :]
sort!(data, :t)
fillmissing!(data)
@assert data == dropmissing(data) # 데이터 무결성 검사
dropmissing!(data)
trng = data[data.t .< Date(2021, 1, 1), :]
test = data[data.t .≥ Date(2020, 12, 31), :]

selected = [:WTI유, :구리, :금]

# scatter(trng[:, :WTI유], trng[:, :천연가스])


# X = Matrix(data[:, Not(:t)])'
X = Matrix(trng[:, selected])'
DX = diff(X, dims = 2); X = X[:, 1:(end-1)]
ddprob = ContinuousDataDrivenProblem(X, DX)

@variables u[1:size(X)[1]]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)

opt = STLSQ(10^(-4), 10)
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 4))
soleq = get_basis(ddsol);
get_parameter_map(soleq)
print(soleq, true) # soleq[1].rhs
is_converged(ddsol)

f̂ = dynamics(soleq)
P = get_parameter_values(soleq)

@time using OrdinaryDiffEq
u0 = X[:, end]
tspan = (30, 50)
dt = 1
DDM = solve(ODEProblem(f̂, u0, tspan, P), RK4(), saveat = dt)

@time using Plots
p1 = plot(DDM, idxs = (0,1))
p2 = plot(DDM, idxs = (0,2))
p3 = plot(DDM, idxs = (0,3))
plot!(p1, 0:30, X[1, (end-30):end], xlims = (0,50), color = :black, label = :none, title = "WTI")
plot!(p2, 0:30, X[2, (end-30):end], xlims = (0,50), color = :black, label = :none, title = "Copper")
plot!(p3, 0:30, X[3, (end-30):end], xlims = (0,50), color = :black, label = :none, title = "Gold")
plot!(p1, 30:50, test[1:21, 2], color = :black, label = "Data")
plot!(p2, 30:50, test[1:21, 3], color = :black, label = "Data")
plot!(p3, 30:50, test[1:21, 4], color = :black, label = "Data")
plot(p1, p2, p3, layout = (3, 1), size = (800, 600))

### ---

@time using DataDrivenDiffEq
@time using ModelingToolkit
@time using OrdinaryDiffEq
@time using DataDrivenSparse
@time using LinearAlgebra
@time using Plots

# Create a test problem
function f(vec, p, t)
    x, y, z, u, v, w = vec

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z

    u̇ = - v - w
    v̇ = u + 0.1 * v
    ẇ = 0.1 + w * (u - 14)

    return [ẋ, ẏ, ż, u̇, v̇, ẇ]
end

u0 = [1.0; 0.0; 0.0; 1.0; 1.0; 1.0]
tspan = (0.0, 100.0)
dt = 0.01
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, RK4(), saveat = dt)
plot(
    plot(sol, idxs = (1, 2, 3), color = :black),
    plot(sol, idxs = (4, 5, 6), color = :black)
)
@assert sol(sol.t,Val{1})[1] == ((sol(dt) - sol(0)) / dt)

ddprob = DataDrivenProblem(sol)

vec = @variables x y z u v w
basis = Basis(polynomial_basis(vec, 2), vec)

opt = STLSQ(exp10.(-5:0.1:-1))
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 8))
soleq = get_basis(ddsol)
print(soleq, true) # soleq[1].rhs
get_parameter_map(soleq) # Equlvalent to get_parameter_values(soleq)
ddsol.residuals # Equlvalent to rss(ddsol)


f̂ = dynamics(soleq)
P = get_parameter_values(soleq)
# f̂(u0, P, 0)

DDM = solve(ODEProblem(f̂, u0, tspan, P), RK4(), saveat = dt)
plot(
    plot(sol, idxs = (1, 2, 3), color = :black),
    plot(sol, idxs = (4, 5, 6), color = :black, zlabel = "SimulationModel"),
    plot(DDM, idxs = (1, 2, 3), color = :blue),
    plot(DDM, idxs = (4, 5, 6), color = :blue, zlabel = "DataDrivenModel"),
    size = (600, 600)
)

### ---

@time using OrdinaryDiffEq
@time using DataDrivenSparse
@time using DataDrivenDiffEq
@time using Plots

function f(vec, p, t)
    x, y, z = vec

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z

    return [ẋ, ẏ, ż]
end
p,t = Nothing, Nothing
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
dt = 0.01
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, RK4(), saveat = dt)

X = Array(sol)
DX1 = hcat(sol.prob.f.(eachcol(X), p, t)...)
DX2 = Array(sol(sol.t, Val{1}))

ddprob1 = ContinuousDataDrivenProblem(X, DX1)
ddprob2 = ContinuousDataDrivenProblem(X, DX2)

@variables u[1:size(X)[1]]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)
opt = STLSQ(exp10.(-5:0.1:-1), 0)

ddsol1 = solve(ddprob1, basis, opt, options = DataDrivenCommonOptions(digits = 6))
ddsol2 = solve(ddprob2, basis, opt, options = DataDrivenCommonOptions(digits = 6))

soleq1 = get_basis(ddsol1); get_parameter_values(soleq1)
soleq2 = get_basis(ddsol2); get_parameter_values(soleq2)

f̂1 = dynamics(soleq1);
P1 = get_parameter_values(soleq1)
f̂2 = dynamics(soleq2);
P2 = get_parameter_values(soleq2)

DDMsol1 = solve(ODEProblem(f̂1, u0, tspan, P1), RK4(), saveat = dt);
DDMsol2 = solve(ODEProblem(f̂2, u0, tspan, P2), RK4(), saveat = dt);

plot(
    plot(sol, idxs = (1,2,3), title = "Ground Truth"),
    plot(DDMsol1, idxs = (1,2,3), title = "Exact"),
    plot(DDMsol2, idxs = (1,2,3), title = "Order 1"),
    size = (600, 600)
)

rss(ddsol1)
rss(ddsol2)

using FiniteDifferences

central_fdm(5, 1)
backward_fdm(2, 1)

include("../src/Utils.jl")
fdiff(X, stencil = 2, dt = dt)
DX2[:, 2:end]

fdiff(X, stencil = 5, dt = dt)
DX1[:, 5:end]
DX2[:, 5:end]

sum(abs, DX1[:, 5:end] - fdiff(X, order = 5, dt = dt))
sum(abs, DX1[:, 5:end] - DX2[:, 5:end])

# https://symbolics.juliasymbolics.org/dev/manual/sparsity_detection/
bart = Dict(get_parameter_map(soleq1)) |> keys
foo = Symbolics.arguments(Symbolics.arguments(soleq1[2].rhs)[3])[1]
foo in bart

### ---

@time using ProgressBars, StatsBase
@time using OrdinaryDiffEq
@time using DataDrivenSparse
@time using DataDrivenDiffEq
@time using Plots
include("../src/Utils.jl")

function f(vec, p, t)
    x, y, z = vec

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z

    return [ẋ, ẏ, ż]
end
p,t = Nothing, Nothing
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 150.0)
dt = 0.01
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, RK4(), saveat = dt)[5001:end]

@variables x, y, z
basis = Basis(polynomial_basis([x, y, z], 2), [x, y, z])

# X = Array(sol)[:, rand(1:10001, 8)]
X = Array(sol)
DX = hcat(sol.prob.f.(eachcol(X), p, t)...)

# mse_ = []
# eD_ = []
# eM_ = []
# for points in ProgressBar(2:2:20)
# # for k in ProgressBar(1:100)
#     rX = X[:, 1:k:end]
#     rDX = DX[:, 1:k:end]
#     # println(size(rX)[2])
#     DY, Y = fdiff(rX, stencil = points, dt = k*dt, method = central_fdm); DY
#     # sum(abs, Y - rX[:, (1+ (points ÷ 2)):(end- (points ÷ 2))])
#     push!(eD_, mean(abs, DY - rDX[:, (1+ (points ÷ 2)):(end- (points ÷ 2))]))

#     ddprob = ContinuousDataDrivenProblem(Y, DY)
#     # ddprob = ContinuousDataDrivenProblem(X, DX)

#     opt = STLSQ(exp.(-3:1:0), 0)
#     # opt = STLSQ(exp10.(-15:0.1:5), 0)
#     ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 6))
#     soleq = get_basis(ddsol)
#     f̂ = dynamics(soleq)
#     P = get_parameter_values(soleq)
    
#     # print(soleq); get_parameter_map(soleq)
#     push!(mse_, rss(ddsol) / length(Y))
#     push!(eM_, sum(abs, f̂([-1,0,-1], P, 0) - [10, -29, 8/3]))
# end
# plot(
#     plot(mse_, ylabel = "MSE"),
#     scatter(eM_, mse_),
#     plot(eD_, xlabel = "stencil points(x2)", ylabel = "Error of derivatives"),
#     plot(eM_, ylabel = "My error"),
#     plot_title = "Lack of data = 100 (fixed)", layout = (2,2), legend = :none
# )

mse__ = zeros(10, 20)
mye__ = zeros(10, 20)
for (i, points) in enumerate(2:2:20)
    for (j, k) in enumerate(1:20) |> ProgressBar
        rX = X[:, 1:k:end]
        rDX = DX[:, 1:k:end]
        # println(size(rX)[2])
        DY, Y = fdiff(rX, stencil = points, dt = k*dt, method = central_fdm); DY
        # sum(abs, Y - rX[:, (1+ (points ÷ 2)):(end- (points ÷ 2))])
        push!(eD_, mean(abs, DY - rDX[:, (1+ (points ÷ 2)):(end- (points ÷ 2))]))
    
        ddprob = ContinuousDataDrivenProblem(Y, DY)
        # ddprob = ContinuousDataDrivenProblem(X, DX)
    
        opt = STLSQ(exp.(-3:1:0), 0)
        # opt = STLSQ(exp10.(-15:0.1:5), 0)
        ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 6))
        soleq = get_basis(ddsol)
        f̂ = dynamics(soleq)
        P = get_parameter_values(soleq)
        
        # print(soleq); get_parameter_map(soleq)
        # push!(mse_, rss(ddsol) / length(Y))
        # push!(eM_, sum(abs, f̂([-1,0,-1], P, 0) - [10, -29, 8/3]))
        mse__[i,j] = rss(ddsol) / length(Y)
        mye__[i,j] = sum(abs, f̂([-1,0,-1], P, 0) - [10, -29, 8/3])
    end
end
Plots.plotly()
plot(1:20, 1:10, mse__)
plot(1:20, 1:10, mse__, st = :surface, title = "dt = 0.01", xlabel = "Lack of data", ylabel = "stencil points(x2)", zlabel = "MSE", size = (800, 600)); savefig("dt = 0.01, MSE.html")
plot(1:20, 1:10, mye__, st = :surface, title = "dt = 0.01", xlabel = "Lack of data", ylabel = "stencil points(x2)", zlabel = "My Error", size = (800, 600)); savefig("dt = 0.01, My Error.html")

Plots.gr()
heatmap(mse__, size = (600, 600), title = "MSE, dt = 0.01", xlabel = "Lack of data", ylabel = "stencil points(x2)", zlabel = "MSE")
heatmap(mye__, size = (600, 600), title = "My error, dt = 0.01", xlabel = "Lack of data", ylabel = "stencil points(x2)", zlabel = "My Error")

png("temp")

### ---

@time using CSV, DataFrames, Dates

function fillmissing!(data)
    for col in eachcol(data)
        while true
            bit_missing = ismissing.(col)
            if sum(bit_missing) == 0 break end
            col[bit_missing] .= ((circshift(col, -1) + circshift(col, 1)) / 2)[bit_missing]
            bit_missing = ismissing.(col)
            col[bit_missing] .= circshift(col, 1)[bit_missing]
        end
    end
end

function mydata()
    data = DataFrame()
    DATA_ = []
    for fn = readdir("data")
        data
        itemname = Symbol(first(split(fn, " ")))
        if (itemname == :알루미늄) || (itemname == :오렌지) continue end
        push!(DATA_, CSV.read("data/" * fn, DataFrame))
        DATA_[end].날짜 = Date.(DATA_[end].날짜, dateformat"y- m- d")
        select!(DATA_[end], ["날짜", "종가"])
        rename!(DATA_[end], [:t, itemname])
        if eltype(DATA_[end][:, 2]) <: AbstractString
            DATA_[end][:, 2] = replace.(DATA_[end][:, 2], "," => "")
            DATA_[end][!, 2] = parse.(Float64, DATA_[end][:, 2])
        end
        if isempty(data)
            data = deepcopy(DATA_[end])
        else
            leftjoin!(data, DATA_[end], on = :t, makeunique = true)
        end
    end
    data = data[data.t .< Date(2022, 1, 1), :]
    sort!(data, :t)
    fillmissing!(data)
    @assert data == dropmissing(data) # 데이터 무결성 검사
    dropmissing!(data)
    return data
end

data = mydata()
data = data[data.WTI유 .> 0, :]

for col in eachcol(data)[Not(1)]
    col .= (col ./ maximum(col)) #* 2π
end

trng = data[data.t .< Date(2021, 1, 1), :]
test = data[data.t .≥ Date(2020, 12, 31), :]

@time using OrdinaryDiffEq
@time using DataDrivenDiffEq
@time using DataDrivenSparse
@time using Plots; default(fontfamily = "")

selected = [:WTI유, :구리]
X = Matrix(trng[:, Not(:t)])'
# X = Matrix(trng[:, selected])'
# DX = diff(X, dims = 2); X = X[:, 1:(end-1)]; DX
DX, X = fdiff(X, stencil = 2, dt = 1, method = central_fdm)
m, n = size(X)

@time using NoiseRobustDifferentiation
tvDX = hcat([tvdiff(x, 100, 1000, dx = 1) for x in eachrow(X)]...)
plot(
    plot(X[1,:], title = "Data(WTI)"),
    plot(tvDX'[1,:], title = "TVdiff"),
    plot(DX[1,:], title = "FDM"), layout = (3,1), size = (1000, 1000), legend = :none
)

# ddprob = ContinuousDataDrivenProblem(X, DX)
ddprob = ContinuousDataDrivenProblem(X, tvDX')

@variables u[1:(size(X)[1])]
u = DataDrivenDiffEq.scalarize(u)
basis = Basis(polynomial_basis(u, 2), u)
# basis = Basis(fourier_basis(u, 20), u)
# basis = Basis(sin_basis([i - j for j in u for i in u], 1), u)
# basis = Basis([fourier_basis(u, 20); polynomial_basis(u, 2)], u)

opt = STLSQ(10^(-6), 10)
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 6))
soleq = get_basis(ddsol);
get_parameter_map(soleq)
print(soleq, true) # soleq[1].rhs
is_converged(ddsol)

f̂ = dynamics(soleq);
P = get_parameter_values(soleq)

u0 = collect(trng[end, Not(1)])
tend = 250; tspan = (30, tend)
dt = 1
DDM = solve(ODEProblem(f̂, u0, tspan, P), Tsit5(), saveat = dt)

plt_ = []
for k in 1:size(X)[1]
    ptemp = plot(DDM, idxs = (0,k), legend = :none)
    k = k+1
    plot!(ptemp, 0:30, trng[(end-30):end, k], xlims = (0,tend), color = :black, title = names(data)[k])
    plot!(ptemp, 30:tend, test[1:(tend-29), k], color = :black)
    push!(plt_, ptemp)
end
plot(plt_..., size = (1000,1000), layout = (:,2))

plot(vec(sum(X, dims = 1)))
plot(vec(sum(DX, dims = 1)))

histogram(vec(sum(DX, dims = 1)))


using Interpolations

y = vec(rand(10))
cubic_f = cubic_spline_interpolation(1:10, y)
plot(1:0.1:10, cubic_f.(1:0.1:10), size = (800, 800)); scatter!(1:10, y)

Y_ = []
for i in 1:m
    y = vec(X[i, :])
    cubic_f = cubic_spline_interpolation(1:4610, y)
    push!(Y_, cubic_f.(1:0.1:4610))
end
Y = hcat(Y_...)'
DY, Y = fdiff(Y, stencil = 20, dt = 0.1, method = central_fdm)
ddprob = ContinuousDataDrivenProblem(Y, DY)

plot(1:0.1:4610, cubic_f.(1:0.1:4610), size = (800, 800))
scatter!(1:4610, y, msw = 0, ma = 0.5, ms = 2)

# function f(vec, p, t)
#     x, y, z, u, v, w = vec

#     ẋ = 10.0 * (y - x)
#     ẏ = x * (28.0 - z) - y
#     ż = x * y - (8 / 3) * z

#     u̇ = - v - w
#     v̇ = u + 0.1 * v
#     ẇ = 0.1 + w * (u - 14)

#     return [ẋ, ẏ, ż, u̇, v̇, ẇ]
# end
function f(vec, p, t)
    x, y, z = vec

    ẋ = 10.0 * (y - x)
    ẏ = x * (28.0 - z) - y
    ż = x * y - (8 / 3) * z

    return [ẋ, ẏ, ż]
end

p,t = Nothing, Nothing
u0 = [1; 0; 0]
u0 = [-8; 8; 27]
tspan = (0.0, 200.0)
dt = 0.001
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, RK4(), saveat = dt)
plot(sol, idxs = (1, 2, 3), color = :black)
@assert sol(sol.t,Val{1})[1] == ((sol(dt) - sol(0)) / dt)

ddprob = DataDrivenProblem(sol)

using FiniteDifferences
include("Utils.jl")
X = Array(sol)
DX, X = fdiff(X, stencil = 2, dt = dt, method = central_fdm)
# DX = hcat(sol.prob.f.(eachcol(X), p, t)...)
ddprob = ContinuousDataDrivenProblem(X, DX)

Y_ = []
for i in 1:3
    y = sol[i,:]
    cubic_f = cubic_spline_interpolation(0:dt:100, y)
    push!(Y_, cubic_f.(0:(0.0001):100))
end
Y = hcat(Y_...)'
DY, Y = fdiff(Y, stencil = 4, dt = 0.0001, method = central_fdm)
ddprob = ContinuousDataDrivenProblem(Y, DY)

xyz = @variables x y z
basis = Basis(polynomial_basis(xyz, 2), xyz)

opt = STLSQ(exp10.(-6:1:-1))
ddsol = solve(ddprob, basis, opt, options = DataDrivenCommonOptions(digits = 6))
soleq = get_basis(ddsol)
print(soleq, true) # soleq[1].rhs
get_parameter_map(soleq) # Equlvalent to get_parameter_values(soleq)
ddsol.residuals # Equlvalent to rss(ddsol)


f̂ = dynamics(soleq)
P = get_parameter_values(soleq)
# f̂(u0, P, 0)

DDM = solve(ODEProblem(f̂, u0, tspan, P), RK4(), saveat = dt)
plot(
    plot(sol, vars = (1, 2, 3), color = :black),
    plot(DDM, vars = (1, 2, 3), color = :blue),
    size = (600, 600)
)

# badpoints = []
for k in 9000:(-10):1
    residuals = [sum(abs2, f̂(X[:,j], P, 0) - DX[:,j]) for j in axes(X)[2]]
    # residuals[badpoints] .= 0
    # push!(badpoints, argmax(residuals))
    # taboo = Not(badpoints)
    taboo = ordinalrank(residuals) .< k
    newddprob = ContinuousDataDrivenProblem(X[:, taboo], DX[:, taboo])
    newddsol = solve(newddprob, basis, opt, options = DataDrivenCommonOptions(digits = 4))
    newsoleq = get_basis(newddsol)
    f̂ = dynamics(newsoleq)
    P = get_parameter_values(newsoleq)
        println(get_parameter_map(newsoleq))
        println("k = ", k, "\trss: ", rss(newddsol))
    # print(newsoleq, true)
end

using StatsBase
x = rand(10)

