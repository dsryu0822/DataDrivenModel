using Base.Threads: @threads # Base.Threads.nthreads()
using DataFrames, CSV, ProgressMeter

function RK4(f::Function, v::AbstractVector, h=10^(-2))
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end

function factory_lorenz(idx::Int64, ρ::Number; ic = [10.,10.,10.], tspan = [0., 10.])
    σ = 10
    β = 8/3
    function lorenz(v::AbstractVector)
        x, y, z = v
        dx = σ*(y - x)
        dy = x*(ρ - z) - y
        dz = x*y - β*z
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


r_range = 0:0.1:100
schedules = DataFrame(idx = eachindex(r_range), r = r_range)
# CSV.write("G:/DDM/lyapunov/lorenz_schedules.csv", schedules)
# schedules = CSV.read("G:/DDM/lyapunov/lorenz_schedules.csv", DataFrame)
@showprogress @threads for dr = eachrow(schedules)
    data = factory_lorenz(DataFrame, dr.idx, dr.r, tspan = [0, 1000])[500001:end, :]
    CSV.write("G:/DDM/lyapunov/_lorenz/$(lpad(dr.idx, 5, '0')).csv", data)
end