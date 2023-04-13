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