using DecisionTree, Random, StatsBase, Dates; @info now()
using Base.Threads: @threads # Base.Threads.nthreads()
include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")

J_lorenz(x,y,z,ρ) = [
     -10   10  0
      ρ-z -1  -x
        y  x  -8/3
]
function lyapunov_lorenz_GS()
    r_range = 0:0.1:100
    schedules = DataFrame(idx = eachindex(r_range), r = r_range)
    result = DataFrame(ρ = Float64[], λ1 = Float64[], λ2 = Float64[], λ3 = Float64[])
    dt = 1e-3
    @showprogress @threads for dr = eachrow(schedules)
        filename = "lyapunov/lorenz/$(lpad(dr.idx, 5, '0')).csv"
        !isfile(filename) && continue
        
        λ = zeros(3);
        data = CSV.read(filename, DataFrame)
        J = J_lorenz(collect(data[1, 1:3])..., dr.r)
        U = I(3)
        for i in 2:nrow(data)
            U, V = gram_schmidt(U)
            λ += V |> eachcol .|> norm .|> log          

            U = RK4(J, U, dt)
            J = J_lorenz(collect(data[i, 1:3])..., dr.r)
        end

        λ ./= dt*nrow(data)
        push!(result, [dr.r, λ...])
    end
    sort!(result, :ρ)
    CSV.write("lyapunov/lorenzGS.csv", result, bom = true)
    plot(xticks = 0:20:100, legend = :none, size = [600, 300])
    plot!(result.ρ, result.λ1, lw = 2, color = 1)
    plot!(result.ρ, result.λ2, lw = 2, color = 2)
    plot!(result.ρ, result.λ3, lw = 2, color = 3)
    png("lyapunovGS")
end
lyapunov_lorenz_GS()

function lyapunov_lorenz_QR()
    r_range = 0:0.1:100
    schedules = DataFrame(idx = eachindex(r_range), r = r_range)
    result = DataFrame(ρ = Float64[], λ1 = Float64[], λ2 = Float64[], λ3 = Float64[])
    dt = 1e-3
    @showprogress @threads for dr = eachrow(schedules)
        filename = "lyapunov/lorenz/$(lpad(dr.idx, 5, '0')).csv"
        !isfile(filename) && continue
        
        λ = zeros(3);
        data = CSV.read(filename, DataFrame)
        J = J_lorenz(collect(data[1, 1:3])..., dr.r)
        U, _ = qr(J); U = Matrix(U)
        for i in 2:nrow(data)
            U, V = gram_schmidt(U)
            λ += V |> eachcol .|> norm .|> log          

            U = RK4(J, U, dt)
            J = J_lorenz(collect(data[i, 1:3])..., dr.r)
        end

        λ ./= dt*nrow(data)
        push!(result, [dr.r, λ...])
    end
    sort!(result, :ρ)
    CSV.write("lyapunov/lorenzQR.csv", result, bom = true)
    plot(xticks = 0:20:100, legend = :none, size = [600, 300])
    plot!(result.ρ, result.λ1, lw = 2, color = 1)
    plot!(result.ρ, result.λ2, lw = 2, color = 2)
    plot!(result.ρ, result.λ3, lw = 2, color = 3)
    png("lyapunovQR")
end
lyapunov_lorenz_QR()