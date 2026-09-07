include.("../core/" .* readdir("core")[[1,2,3,4,6]])

using FFTW

function denoise(x)
    fftz = fft(x)
    # θ = sort(abs.(fftz), rev = true)[8]
    θ = percentile(abs.(fftz), 95)
    tfftz = (θ .< abs.(fftz)) .* fftz
    return abs.(ifft(tfftz))
end

nino34 = CSV.read("G:/seasurface/nino34Tlong.csv", DataFrame)

sort(abs.(fft(nino34.k21[t_[2024]])))

plot(nino34.date, nino34.k4121, size = [800, 200], xformatter = _ -> "")
plot(nino34.date, denoise(nino34.k4121), size = [800, 200], xformatter = _ -> "")
plot(nino34.k4121[1:365])

mSST = mean(Matrix(nino34[:, Not(:date)]), dims = 2)[:]
μ = mean(mSST)
data = add_diff(rename(add_diff(DataFrame(u = (denoise(mSST) .- μ)), method = :TVD), ["y", "x"]), method = :TVD)
plot(
    plot(nino34.date, data.x),
    plot(nino34.date, data.y, yticks = [0]),
    plot(nino34.date, data.dy, yticks = [0]),
    xticks = Date(2002):Year(1):Date(2026),
    size = [1000, 600], layout = (:, 1), xformatter = x -> string(Int64(year(x)))
)

t_ = Dict([t => findall(Date(t) .≤ nino34.date .< Date(t + 1)) for t = 2003:2025])
nino34_ = Dict([t => data.x[Date(t) .≤ nino34.date .< Date(t + 1)] for t = 2003:2025])
nino34d = Dict([t => data.y[Date(t) .≤ nino34.date .< Date(t + 1)] for t = 2003:2025])
plot(nino34_[2022] .+ μ)
plot(
    plot(nino34_[2003]),
    plot(nino34_[2003], nino34d[2003]),
    size = [800, 400]
)


vrbl = (["dy", "dx"], ["y", "x"])

cnfg = cook(vrbl, poly = 0:2)
λ_ = exp10.(-8:-1)
oλ_ = Dict(2003:2025 .=> 0.0)
for t = keys(t_)
    mse_ = fill(Inf, length(λ_))
    for k in eachindex(λ_)
        f0 = SINDy(data[t_[t], :], vrbl, cnfg; λ = λ_[k]);
        traj0 = ssolve(f0, data[t_[t], last(vrbl)][1, :], 0:1e-1:(length(t_[t])))[1:10:(end-1), :]
        if nrow(traj0) == length(t_[t])
            mse_[k] = sum((traj0.x .- data[t_[t], :x]).^2)
        end
    end
    oλ_[t] = λ_[argmin(mse_)]
    f0 = SINDy(data[t_[t], :], vrbl, cnfg; λ = oλ_[t]); f0 |> print
    traj0 = ssolve(f0, data[t_[t], last(vrbl)][1, :], 0:1e-1:(365))[1:10:end, :]
    plot(data[t_[t], :x])
    plot!(traj0.x)
    png("SST_$t")
end


# ✅ 2020, 2025
y0 = 2016
y1 = 2025
colnames = ["y0", "y1", "nfails"]
result = DataFrame([[] for _ in colnames], colnames)
for y0 = 2003:2025
    f0 = SINDy(data[t_[y0], :], vrbl, cnfg; λ = oλ_[y0]); f0 |> print
    for y1 = 2003:2025
        f1 = SINDy(data[t_[y1], :], vrbl, cnfg; λ = oλ_[y1]); f1 |> print
        if y0 ≥ y1 continue end

        plt = plot(xticks = [[2000:5:2030]...; [y0, y1]])
        gt = [maximum(nino34_[t]) + μ for t = 2003:2025]
        # gt = [maximum(data[t_[t], :x]) + μ for t = 2003:2025]
        scatter!(plt, 2003:2025, gt, color = :black, msw = 0)
        plot!(plt, 2003:2025, gt, color = :black, msw = 0)
        vline!(plt, [y0, y1], color = :black)

        # traj0 = ssolve(f0, data[t_[y0], last(vrbl)][1, :], 0:1e-1:365)[1:10:end, :]
        # traj1 = ssolve(f1, data[t_[y1], last(vrbl)][1, :], 0:1e-1:365)[1:10:end, :]

        # plot(nino34_[y0])
        # plot!(traj0.x)
        # plot(nino34_[y1])
        # plot!(traj1.x)

        fβ = affine(Function, f0, f1)

        βm = (2003 - y0) / (y1 - y0)
        βM = (2030 - y0) / (y1 - y0)
        MT = []
        β_ = βm:(1/(y1 - y0)):βM
        y_ = β_*(y1 - y0) .+ y0
        for β = β_
            ic = ((1-β)*[data[t_[y0], last(vrbl)][1, :]...;] + β*[data[t_[y1], last(vrbl)][1, :]...;])
            sol = solve(ODEProblem(fβ, ic, (0, 365)), p = [β], RK4(), dt = 1e-1, adaptive=false, maxiters = Inf)
            push!(MT, maximum(sol[2, :]) + μ)
        end

        bit_failed = isnan.(MT) .|| isinf.(MT)
        push!(result, [y0, y1, sum(bit_failed)])
        if (sum(bit_failed) == 0) && maximum(MT) < 32 && minimum(MT) > 25
        plot(plt, y_, MT, color = :red, ylims = [25, 32], xlims = [2002, 2032], shape = :x)
        png("show $(y0)_$(y1)")
        end
    end
end

# plot(nino34_[2014])
# plot!(prdt2014[2, 1:10:end])

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    denoise

''''''''''''''''''''''''''''''''''''''''''''''''''"""

function denoise(x)
    fftz = fft(x)
    # θ = sort(abs.(fftz), rev = true)[8]
    θ = percentile(abs.(fftz), 99)
    tfftz = (θ .< abs.(fftz)) .* fftz
    return abs.(ifft(tfftz))
end
plot(
    plot(nino34.date, nino34.k4121, xformatter = _ -> "", color = :black),
    plot(nino34.date, denoise(nino34.k4121), xformatter = _ -> "", color = :red),
    layout = (:, 1)
)

plot(
    plot(denoise(nino34.k21[t_[2024]])),
    plot(denoise(nino34.k4121[t_[2024]])),
    plot(denoise(nino34.k8221[t_[2024]])),
    layout = (:, 1)
)
plot(
    plot(denoise(nino34.k21[t_[2024]]), denoise(nino34.k4121[t_[2024]]), denoise(nino34.k8221[t_[2024]])),
    plot(denoise(nino34.k21[t_[2024]]), denoise(nino34.k8221[t_[2024]])),
    size = [800, 400]
)


data_ = Dict([t => add_diff(DataFrame(
    x = denoise(nino34.k21[t_[t]]),
    y = denoise(nino34.k4121[t_[t]]),
    z = denoise(nino34.k8221[t_[t]])),
    method = :TVD) for t in 2003:2025])
vrbl = half(names(data_[2003]))
cnfg = cook(vrbl, poly = 0:2)

for y0 = 2003:2025
    f0 = SINDy(data_[y0], vrbl, cnfg; λ = 1e-16); f0 |> print
    traj0 = ssolve(f0, data_[y0][1, :], 0:1e-1:   365)[1:10:end, :]
    plot(data_[y0].x, data_[y0].y, data_[y0].z)
    plot!(traj0.x, traj0.y, traj0.z)
    png("SST_$y0")
end

y0 = 2016
y1 = 2020
f0 = SINDy(data_[y0], vrbl, cnfg; λ = 1e-2); f0 |> print
f1 = SINDy(data_[y1], vrbl, cnfg; λ = 1e-2); f1 |> print
traj0 = ssolve(f0, data_[y0][1, :], 0:1e-1:365)[1:10:end, :]
traj1 = ssolve(f1, data_[y1][1, :], 0:1e-1:365)[1:10:end, :]

plot(data_[y0].x, data_[y0].y, data_[y0].z)
plot!(traj0.x, traj0.y, traj0.z)
plot(data_[y1].x, data_[y1].y, data_[y1].z)
plot!(traj1.x, traj1.y, traj1.z)

plt = plot(xticks = [[2000:5:2030]...; [y0, y1]])
for t = 2003:2025
    scatter!(plt, [t], [maximum(data_[t].y)], color = :black, msw = 0)
end
vline!(plt, [y0, y1], color = :black)

fβ = affine(Function, f0, f1)

βm = (2003 - y0) / (y1 - y0)
βM = (2030 - y0) / (y1 - y0)
MT = []
β_ = βm:(1/(y1 - y0)):βM
y_ = β_*(y1 - y0) .+ y0
for β = β_
    ic = ((1-β)*[data_[y0][1, last(vrbl)]...] + β*[data_[y1][1, last(vrbl)]...])
    sol = solve(ODEProblem(fβ, ic, (0, 365)), p = [β], RK4(), dt = 1e-1, adaptive=false, maxiters = Inf)
    push!(MT, maximum(sol[2, :]))
end

scatter(plt, y_, MT, color = :red, ylims = [25, 32], xlims = [2002, 2032], shape = :x)

[y_ MT]