include("../core/header.jl")

function factory_lorenz_tipping(ρ0::Number; ic = [0, ρ0, 1.,1.,1.], tspan = [0., 10.], dt = 1e-3)
    σ = 10
    β = 8/3
    ε = 1e-2
    function sys(v::AbstractVector)
        T, ρ, x, y, z = v

        dT = 1
        dρ = -ε
        dx = σ*(y - x)
        dy = x*(ρ - z) - y
        dz = x*y - β*z
        return [dT, dρ, dx, dy, dz]
    end
        
    t_ = first(tspan):dt:last(tspan)
    len_t_ = length(t_)
    
    t, tk = .0, 0
    v = ic; DIM = length(v)
    traj = zeros(len_t_+2, 2DIM)
    while tk ≤ len_t_
        _,ρ,x,y,z = v
        v, dv = RK4(sys, v, dt)

        if t ≥ first(t_)
            tk += 1
            traj[tk+1,         1:DIM ] =  v
            traj[tk  , DIM .+ (1:DIM)] = dv
        end
    end
    return traj[2:(end-2), :]
end
factory_lorenz_tipping(T::Type, args...; kargs...) =
DataFrame(factory_lorenz_tipping(args...; kargs...), ["t", "ρ", "x", "y", "z", "dt", "dρ", "dx", "dy", "dz"])


@time data = factory_lorenz_tipping(DataFrame, 230, tspan = [0., 3000])
# anime = @animate for t0 = 1:10000:(nrow(data)-10000)
#     dat = data[t0 .+ (0:10000), :]
#     plot(dat.x, dat.y, dat.z, xlim = extrema(data.x), ylim = extrema(data.y), zlim = extrema(data.z), legend = :none,
#         title = "t = $(round(dat.t[1])) ~ $(round(dat.t[end]))", color = :black)
# end
# mp4(anime, "temp.mp4")
# plot(data.t[1:100:end], data.ρ[1:100:end])

# scatter(data.t[1:100:end], data.x[1:100:end], color = :black, legend = :none, ms = 1)
vrbl = (names(data)[8:10], names(data)[3:5])
cnfg = cook(last(vrbl), poly = 0:2)

bit_xm = arglmax(data.z)
scatter(data.t[bit_xm], data.z[bit_xm], msw = 0, ylims = [220, 320], ms = 1, color = :black, legend = :none, xlims = [0, 3000])

# static = factory_lorenz(DataFrame, 230)[50000:end, :]
# plot(static.x, static.y, static.z)
# SINDy(static, vrbl, cnfg, λ = 0.01) |> print

hrzn = []
vrtl = []
for ρ in 230:(-0.1):200
    traj = factory_lorenz(DataFrame, ρ, tspan = [0., 50.], dt = 1e-3)[40000:end, :]
    points = arglmax(traj.z)
    push!(vrtl, traj.z[points])
    push!(hrzn, repeat([ρ], length(points)))
end
scatter([hrzn...;], [vrtl...;], xflip = true, xlims = [200, 230], ylims = [220, 320], msw = 0, color = :black, ms = 1, legend = :none, alpha = 0.5)
png("3")

# τ = data.t
# traj0 = data[ 490 .≤ τ .≤  500, :]
# traj1 = data[ 990 .≤ τ .≤ 1000, :]
traj0 = factory_lorenz(DataFrame, 230, tspan = [0., 50.], dt = 1e-3)[40000:end, :]
traj1 = factory_lorenz(DataFrame, 225, tspan = [0., 50.], dt = 1e-3)[40000:end, :]
plot(
    plot(traj0.x, traj0.y, traj0.z, ticks = [], color = :black),
    plot(traj1.x, traj1.y, traj1.z, ticks = [], color = :black),
    legend = :none
)

@time f0 = SINDy(traj0, vrbl, cnfg, λ = 0.01); f0 |> print
@time f1 = SINDy(traj1, vrbl, cnfg, λ = 0.01); f1 |> print
@info "SINDy done!"
α_ = 0:1e-2:6
traj_ = [[1 1 1.]]
for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg)
    traj = solve(g, last(eachrow(last(traj_))), 0:1e-3:5)[2:end, :]
    push!(traj_, traj)
end
traj_ = vcat(traj_...)
traj_ = [collect(1e-3:1e-3:(1e-3 * size(traj_, 1))) traj_]
rcvd = DataFrame(traj_, [:t, :x, :y, :z])

α_ = 0:1e-2:6
hrzn = []
vrtl = []
for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg)
    traj = solve(g, [1, 1, 1], 0:1e-3:50)[40000:end, :]
    points = arglmax(traj[:, 3])
    push!(vrtl, traj[points, 3])
    push!(hrzn, repeat([α], length(points)))
end
scatter([hrzn...;], [vrtl...;], xlims = [0, 6], ylims = [220, 320], msw = 0, color = :blue, ms = 1, legend = :none, alpha = 0.5); png("3")


plot(rcvd.z[1:100:end])
bit_xm = arglmax(rcvd.z)
scatter(rcvd.t[bit_xm], rcvd.z[bit_xm], msw = 0, ylims = [220, 320], ms = 1, color = :blue, legend = :none, xlims = [0, 3000])

bit_xm = arglmax(rcvd.z)
scatter(rcvd.t[bit_xm], rcvd.z[bit_xm], msw = 0, ylims = [220, 320], ms = 1, color = :blue, legend = :none, xlims = [0, 3000])

png("temp")
