const κ = 400.0
const μ = 172.363

using ProgressBars
packages = [:DataFrames, :CSV, :Plots, :LaTeXStrings]
@time for package in ProgressBar(packages)
    @eval using $(package)
end
println(join(packages, ", "), " loaded!")

include("../src/ODEdata.jl")

function factory_soft(idx::Int64, d::Number)
    d2 = d/2
    
    function soft(tuv)
        t, u, v    = tuv

        impact = ifelse(abs(u) < d2, 0, -(κ^2)*sign(u)*(abs(u)-d2) - μ*v)
        ṫ = 1
        u̇ = v
        v̇ = cospi(t) + impact
        return [ṫ, u̇, v̇]
    end
    dt = 10^(-5); tend = 50
    t_ = 0:dt:tend

    ndatapoints = round(Int64, tend/(10dt))

    len_t_ = length(t_)
    traj = zeros(6, len_t_+1)
    x = [0.0, 0.05853, 0.47898]
    dx = soft(x)
    traj[1:3, 1] = x

    
    for t in 1:length(t_)
        x, dx = RK4(soft, x, dt)
        if t ≥ ndatapoints
            traj[4:6,   t] = dx
            traj[1:3, t+1] = x
        end
    end
    traj = traj[:, (end-ndatapoints):(end-1)]'

    data = DataFrame(traj,
        ["t", "u", "v", "dt", "du", "dv"])

    return data
end

d_range = 0.1:0.0001:0.3
schedule = DataFrame(idx = eachindex(d_range), d = d_range)

xdots = Float64[]; ydots = Float64[]
for dr in ProgressBar(eachrow(schedule))
    data = factory_soft(dr.idx, dr.d)
    
    idx = [false; diff(abs.(data.u) .> (dr.d/2)) .> 0]

    sampled = data.v[idx]
    append!(xdots, fill(dr.d, length(sampled)))
    append!(ydots, sampled)
end
@time a1 = scatter(xdots, ydots,
xlabel = L"d", ylabel = "Impact velocity",
label = :none, msw = 0, color = :black, ms = 0.5, alpha = 0.5, size = (700, 300))

png(a1, "soft_bifurcation")

dr = eachrow(schedule)[501]
data = factory_soft(dr.idx, dr.d)
_data = data[1:1000:end, :]
points = ([data.u data.v data.dv]')

dbscan(points, 0.1)