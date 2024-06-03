include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

# function lyapunov_soft()
#     schedules = CSV.read("G:/DDM/bifurcation/soft_schedules.csv", DataFrame)
#     vrbl = [:dt, :du, :dv], [:t, :u, :v]
#     cnfg = (; f_ = [cospi, sign], λ = 1e-2)
#     dt = 1e-5; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

#     result = DataFrame(d = Float64[], λ1 = Float64[], λ2 = Float64[], λ3 = Float64[])
#     @showprogress @threads for dr = eachrow(schedules)[1:end]
#         try
#             filename = "G:/DDM/bifurcation/soft/$(lpad(dr.idx, 5, '0')).csv"
#             !isfile(filename) && continue
            
#             data = CSV.read(filename, DataFrame)
#             add_subsystem!(data, vrbl, cnfg; (; θ1, θ2, θ3, min_rank)...); # 15~
#             f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(data, :subsystem)]
#             # Dtree = dryad(data, last(vrbl)); # print_tree(Dtree)
#             J_ = jacobian.(f_)

#             λ = zeros(3);

#             # Q = I(3); v = ones(3);
#             J = substitute(J_[data.subsystem[1]], Dict(t => data.t[1]))
#             Q, _ = qr(J)
#             # for _t in 1:dt:tend
#             for i in 2:nrow(data)
#                 # s = apply_tree(Dtree, [_t; v...])
#                 # v, _ = RK4(f_[s], v, dt)
#                 # J = substitute(J_[s], Dict(t => _t))
#                 J = substitute(J_[data.subsystem[i]], Dict(t => data.t[i]))
#                 # Q = (I(3) + dt*J)*Q
#                 Q = Matrix(Q)
#                 V1 = J*Q
#                 V2 = J*(Q + (dt/2)*V1)
#                 V3 = J*(Q + (dt/2)*V2)
#                 V4 = J*(Q + dt*V3)
#                 Q += (dt/6)*(V1 + 2V2 + 2V3 + V4)

#                 Q, R = qr(Q)
#                 λ += R |> diag .|> abs .|> log
#             end

#             λ ./= (tend - 1)
#             # λ ./= dt*nrow(data)
#             push!(result, [dr.d, λ...])
#         catch
#             @warn "Error at $(dr.idx)"
#         end
#     end
#     sort!(result, :d)
#     CSV.write("G:/DDM/lyapunov/soft.csv", result)
#     plot(legend = :none, size = [600, 300])
#     plot!(result.d, result.λ1, lw = 2, color = 1)
#     plot!(result.d, result.λ2, lw = 2, color = 2)
#     plot!(result.d, result.λ3, lw = 2, color = 3)
#     png("G:/DDM/lyapunov/soft_lyapunov.png")
# end
# lyapunov_soft()

function lorenz(v::AbstractVector; ρ = 28)
    x, y, z = v
    dx = 10*(y - x)
    dy = x*(ρ - z) - y
    dz = x*y - (8/3)*z
    return [dx, dy, dz]
end
J_lorenz(x,y,z,ρ) = [
     -10   10  0
      ρ-z -1  -x
        y  x  -8/3
]

function lyapunov_lorenz()
    r_range = 0:0.1:100
    schedules = DataFrame(idx = eachindex(r_range), r = r_range)
    result = DataFrame(ρ = Float64[], λ1 = Float64[], λ2 = Float64[], λ3 = Float64[])
    dt = 1e-3 #; tend = 1_000_000
    @showprogress @threads for dr = eachrow(schedules)[1:1:end]
        filename = "G:/DDM/lyapunov/_lorenz/$(lpad(dr.idx, 5, '0')).csv"
        !isfile(filename) && continue
        
        λ = zeros(3);
        data = CSV.read(filename, DataFrame)
        J = J_lorenz(collect(data[1, 1:3])..., dr.r)
        Q, _ = qr(J)
        for i in 2:nrow(data)
            J = J_lorenz(collect(data[i, 1:3])..., dr.r)

            Q = Matrix(Q)
            V1 = J*Q
            V2 = J*(Q + (dt/2)*V1)
            V3 = J*(Q + (dt/2)*V2)
            V4 = J*(Q + dt*V3)
            Q += (dt/6)*(V1 + 2V2 + 2V3 + V4)

            Q, R = qr(Q)
            λ += R .|> abs |> diag .|> log
        end

        λ ./= dt*nrow(data)
        push!(result, [dr.r, λ...])
    end
    sort!(result, :ρ)
    CSV.write("G:/DDM/lyapunov/lorenz.csv", result)
    plot(xticks = 0:20:100, legend = :none, size = [600, 300])
    plot!(result.ρ, result.λ1, lw = 2, color = 1)
    plot!(result.ρ, result.λ2, lw = 2, color = 2)
    plot!(result.ρ, result.λ3, lw = 2, color = 3)
    png("lyapunovQR")
end
lyapunov_lorenz()