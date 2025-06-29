include("../core/header.jl")

vrbl = [:dx, :dy, :dz], [:x, :y, :z]
cnfg = (; N = 2)
start = sort(rand(1:nrow(data), 10))

datasets28 = [factory_lorenz(DataFrame, 28, tspan = [0, 100])[(1:100:10000) .+ stt, :] for stt in start]

plot(data.x[1:10:end])

D = zeros(length(datasets28), length(datasets28))
for (i,j) in combinations(1:length(datasets28), 2)
    if i > j continue end
    f = SINDy([datasets28[i]; datasets28[j]], vrbl...; cnfg...)
    D[i,j] = f.MSE
    D[j,i] = f.MSE
end

datasets29 = [factory_lorenz(DataFrame, 29, tspan = [0, 100])[(1:100:10000) .+ stt, :] for stt in start]
datasets40 = [factory_lorenz(DataFrame, 40, tspan = [0, 100])[(1:100:10000) .+ stt, :] for stt in start]

datasets = [datasets28; datasets29; datasets40]
D = fill(Inf, length(datasets), length(datasets))
for (i,j) in combinations(1:length(datasets), 2)
    if i > j continue end
    f = SINDy([datasets[i]; datasets[j]], vrbl...; cnfg...)
    D[i,j] = f.MSE
    D[j,i] = f.MSE
end
M11 = fill(:red, 10, 10)
M22 = fill(:blue, 10, 10)
M33 = fill(:green, 10, 10)
M00 = fill(:black, 10, 10)
colorM = [M11 M00 M00;
          M00 M22 M00;
          M00 M00 M33]
scatter(vec(D), yscale = :log10, legend = :none, color = vec(colorM), msw = 0)
