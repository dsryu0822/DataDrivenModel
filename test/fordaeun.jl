
include("../core/header.jl")

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
data = factory_lorenz(DataFrame, 0, 28, tspan = [0, 100])[50000:end, :]

SINDy(data, [:dx, :dy, :dz], [:x, :y, :z], N = 2)