include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

using LinearAlgebra, StatsBase

# ρ = 28
function lorenz(v::AbstractVector; ρ = 28)
    x, y, z = v
    dx = 10*(y - x)
    dy = x*(ρ - z) - y
    dz = x*y - 8/3*z
    return [dx, dy, dz]
end
J_lorenz(x,y,z,ρ) = [
     -10 10  0
    -z+ρ -1 -x
       y  x -8/3
]

function factory_lorenz(idx::Int64, ρ::Number; ic = [10.,10.,10.], tspan = [0., 10.])
    ___σ = 10
    ___β = 8/3
    function lorenz(v::AbstractVector)
        x, y, z = v
        dx = ___σ*(y - x)
        dy = x*(ρ - z) - y
        dz = x*y - ___β*z
        return [dx, dy, dz]
    end

    dt = 10^(-3)
    t_ = first(tspan):dt:last(tspan)
    
    ndatapoints = count(first(tspan) .< t_ .≤ last(tspan))
    len_t_ = length(t_)

    v = ic; DIM = length(v)
    traj = zeros(2DIM, len_t_+1)
    traj[1:DIM, 1] = v

    for tk in eachindex(t_)
        v, dv = RK4(lorenz, v, dt)
        if tk+1 ≥ (len_t_ - ndatapoints)
            traj[        1:DIM , tk+1] =  v
            traj[DIM .+ (1:DIM), tk  ] = dv
        end
    end
    traj = traj[:, 1:(end-1)]'
    traj = traj[(end-ndatapoints):end, :]

    return traj
end
factory_lorenz(T::Type, args...; ic = [10,10,10.], tspan = [0., 100.]) =
DataFrame(factory_lorenz(args...; ic = ic, tspan), ["x", "y", "z", "dx", "dy", "dz"])

# r_range = 25:0.1:325
# schedules = DataFrame(idx = eachindex(r_range), r = r_range)
# CSV.write("G:/DDM/lyapunov/lorenz_schedules.csv", schedules)
schedules = CSV.read("G:/DDM/lyapunov/lorenz_schedules.csv", DataFrame)
# dr = eachrow(schedules)[1]
@showprogress @threads for dr in eachrow(schedules)
    data = factory_lorenz(DataFrame, dr.idx, dr.r, tspan = [0, 100])[90000:end, :]
    CSV.write("G:/DDM/lyapunov/lorenz/$(lpad(dr.idx, 5, '0')).csv", data)
end

# plot(data.z)
# sampledz = findall((circshift(data.z, 1) .< data.z) .&& (circshift(data.z, -1) .< data.z))
# scatter!(    sampledz
# , data.z[sampledz]
# )

idcs = Int64[]; vrtc = Float64[]; hrzn = Float64[]
@showprogress for dr = eachrow(schedules)
    filename = "G:/DDM/lyapunov/lorenz/$(lpad(dr.idx, 5, '0')).csv"
    !isfile(filename) && continue

    data = CSV.read(filename, DataFrame)
    idx_sampled = findall((circshift(data.z, 1) .< data.z) .&& (circshift(data.z, -1) .< data.z))[2:(end-1)]
    # pop!(idx_sampled)
    sampledz = data.z[idx_sampled]
    append!(idcs, fill(dr.idx, length(sampledz)))
    append!(hrzn, fill(dr.r, length(sampledz)))
    append!(vrtc, sampledz)
end
CSV.write("G:/DDM/lyapunov/lorenz_bifurcation.csv", DataFrame(; vrtc, hrzn))
scatter(hrzn, vrtc, color = :black, ms = 1, legend = :none, msw = 0, ma = 0.1, xticks = [25, 325], size = [800, 800], xlabel = L"\rho");
png("G:/DDM/lyapunov/lorenz_bifurcation.png")

r_range = 0:0.1:100
schedules = DataFrame(idx = eachindex(r_range), r = r_range)
hrzn = Float64[]; vrtc1, vrtc2, vrtc3 = Float64[], Float64[], Float64[]
h = 0.001; tend = 100000
@showprogress for dr in eachrow(schedules)[1:10:end]
    # filename = "G:/DDM/lyapunov/lorenz/$(lpad(dr.idx, 5, '0')).csv"
    # !isfile(filename) && continue
    # data = CSV.read(filename, DataFrame)

    U = I(3); v = ones(3); λ = zeros(3);
    for i in 1:tend
        J = J_lorenz(v..., dr.r)
        U = (I(3) + h*J)*U
        U, R = qr(U)
        v = v + h*lorenz(v, ρ = dr.r)
        λ += R .|> abs |> diag .|> log
    end
    push!(hrzn, dr.r); push!(vrtc1, λ[1] / tend); push!(vrtc2, λ[2] / tend); push!(vrtc3, λ[3] / tend)
end
plot(xticks = 0:20:100, legend = :none)
plot!(hrzn, vrtc1/h, colozr = 1)
plot!(hrzn, vrtc2/h, color = 2)
plot!(hrzn, vrtc3/h, color = 3)