include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

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

# function lorenz(v::AbstractVector; ρ = 28)
#     x, y, z = v
#     dx = 10*(y - x)
#     dy = x*(ρ - z) - y
#     dz = x*y - (8/3)*z
#     return [dx, dy, dz]
# end
# J_lorenz(x,y,z,ρ) = [
#      -10   10  0
#       ρ-z -1  -x
#         y  x  -8/3
# ]

# function lyapunov_lorenz()
#     r_range = 0:0.1:100
#     schedules = DataFrame(idx = eachindex(r_range), r = r_range)
#     result = DataFrame(ρ = Float64[], λ1 = Float64[], λ2 = Float64[], λ3 = Float64[])
#     h = 1e-3 #; tend = 1_000_000
#     @showprogress @threads for dr = eachrow(schedules)[1:1:end]
#         filename = "G:/DDM/lyapunov/lorenz/$(lpad(dr.idx, 5, '0')).csv"
#         !isfile(filename) && continue
        
#         λ = zeros(3);
#         data = CSV.read(filename, DataFrame)
#         J = J_lorenz(collect(data[1, 1:3])..., dr.r)
#         Q, _ = qr(J)
#         for i in 2:nrow(data)
#             J = J_lorenz(collect(data[i, 1:3])..., dr.r)
#             Q = (I(3) + h*J)*Q
#             Q, R = qr(Q)
#             λ += R .|> abs |> diag .|> log
#         end

#         λ ./= h*nrow(data)
#         push!(result, [dr.r, λ...])
#     end
#     sort!(result, :ρ)
#     CSV.write("G:/DDM/lyapunov/lorenz.csv", result)
#     plot(xticks = 0:20:100, legend = :none, size = [600, 300])
#     plot!(result.ρ, result.λ1, lw = 2, color = 1)
#     plot!(result.ρ, result.λ2, lw = 2, color = 2)
#     plot!(result.ρ, result.λ3, lw = 2, color = 3)
#     png("lyapunov")
# end
# lyapunov_lorenz()

# vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi, sign], λ = 1e-2)
# # cnfg = (; M = 3, λ = 1e-2)
# dt = 1e-5; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

# data = CSV.read("G:/DDM/bifurcation/soft/00001.csv", DataFrame)

# add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 15~
# f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
# Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)

# print(f_[1])
# print(f_[2])
# jacobian.(f_)[1]
# jacobian.(f_)[2]


function lyapunov_soft()
    schedules = CSV.read("G:/DDM/bifurcation/soft_schedules.csv", DataFrame)
    vrbl = [:dt, :du, :dv], [:t, :u, :v]
    cnfg = (; f_ = [cospi, sign], λ = 1e-2)
    dt = 1e-5; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
    
    result = DataFrame(d = Float64[], λ1 = Float64[], λ2 = Float64[], λ3 = Float64[])
    @showprogress @threads for dr = eachrow(schedules)[1:end]
        try
            filename = "G:/DDM/bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv"
            !isfile(filename) && continue
            
            data = CSV.read(filename, DataFrame)
            add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 15~
            f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
            # Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
            J_ = jacobian.(f_)

            λ = zeros(3);
            J = substitute(J_[data.subsystem[1]], Dict(t => data.t[1]))
            Q, _ = qr(J)
            for i in 1:nrow(data)
                J = substitute(J_[data.subsystem[i]], Dict(t => data.t[i]))
                Q = (I(3) + dt*J)*Q
                Q, R = qr(Q)
                λ += R .|> abs |> diag .|> log
            end

            λ ./= dt*nrow(data)
            push!(result, [dr.d, λ...])
        catch
            @warn "Error at $(dr.idx)"
        end
    end
    sort!(result, :d)
    CSV.write("G:/DDM/lyapunov/soft.csv", result)
    plot(legend = :none, size = [600, 300])
    plot!(result.d, result.λ1, lw = 2, color = 1)
    plot!(result.d, result.λ2, lw = 2, color = 2)
    plot!(result.d, result.λ3, lw = 2, color = 3)
    png("G:/DDM/lyapunov/soft_lyapunov.png")
end
lyapunov_soft()