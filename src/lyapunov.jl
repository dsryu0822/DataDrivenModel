using DecisionTree, Random, StatsBase, Dates; @info now()
using Base.Threads: @threads # Base.Threads.nthreads()
include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")

# J_soft(t, u, v, d) = [
#                0                                0                                 0
#                0                                0                                 1
#      -π*sinpi(t) ifelse(abs(u) ≥ d/2, -160000, 0) ifelse(abs(u) ≥ d/2, -172.363, 0)
# ]
function lyapunov_soft()
    function J_(t, u, v, d)
        _bit = abs(u) ≥ d/2
        __κ2 = 160000
        __μ = 172.363
        return [          0          0         0
                          0          0         1
                -π*sinpi(t) -__κ2*_bit -__μ*_bit]
    end    
    schedules = CSV.read("G:/DDM/lyapunov/soft_schedules_cache.csv", DataFrame)
    # vrbl = [:dt, :du, :dv], [:t, :u, :v]
    # cnfg = (; f_ = [cospi, sign], λ = 1e-2)
    dt = 1e-6; θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;

    # @showprogress @threads for dr = eachrow(schedules)[1:1:end]
    # filename = "G:/DDM/lyapunov/soft/$(lpad(dr.idx, 5, '0')).csv"
    # !isfile(filename) && continue
    # data = CSV.read(filename, DataFrame)
    result = DataFrame(d = Float64[], λ1 = Float64[], λ2 = Float64[], λ3 = Float64[])
    # @showprogress @threads for dr = eachrow(schedules)[883:3:1167]
    @showprogress @threads for dr = eachrow(schedules)[1:1:end]
        data = factory_soft(DataFrame, dr.idx, dr.d; ic = [dr.t, dr.u, dr.v], tspan = [0, 100], dt)[:, 1:3]
        
        λ = zeros(3);
        J = J_soft(collect(data[1, 1:3])..., dr.d)
        U, _ = qr(J); U = Matrix(U)
        for i in 2:nrow(data)
            U, V = gram_schmidt(U)
            λ += V |> eachcol .|> norm .|> log

            U = RK4(J, U, dt)
            J = J_soft(collect(data[i, 1:3])..., dr.d)
        end

        λ ./= dt*nrow(data)
        push!(result, [dr.d, λ...])
    end
    sort!(result, :d)
    CSV.write("G:/DDM/lyapunov/soft_lyapunov.csv", result, bom = true)
    plot(legend = :none, size = [600, 300])
    plot!(result.d, result.λ1, lw = 2, color = 1)
    plot!(result.d, result.λ2, lw = 2, color = 2)
    plot!(result.d, result.λ3, lw = 2, color = 3)
    png("G:/DDM/lyapunov/soft_lyapunov.png")
end
lyapunov_soft()

# plot(data.u[1:100:end])

# bfcn = CSV.read("G:/DDM/bifurcation/soft_bifurcation.csv", DataFrame)
# plt_bfcn = scatter(bfcn.vrtc, bfcn.hrzn, ms = .5, legend = :none);

# lpnv = CSV.read("G:/DDM/lyapunov/soft.csv", DataFrame)
# plt_lpnv = plot(legend = :none);
# plot!(lpnv.d, lpnv.λ1, lw = 2, color = 1);
# plot!(lpnv.d, lpnv.λ2, lw = 2, color = 2);
# plot!(lpnv.d, lpnv.λ3, lw = 2, color = 3);

# plot(plt_bfcn, plt_lpnv, layout = (2, 1), xlims = [0.1881, .2167])
# png("temp")