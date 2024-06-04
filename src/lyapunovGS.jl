include("../core/DDM.jl")
include("../core/factorio.jl")
include("../core/visual.jl")
using DecisionTree, Random, StatsBase
using Base.Threads: @threads # Base.Threads.nthreads()

function gram_schmidt(J)
    N = size(J, 1)
    U, V = deepcopy(J), deepcopy(J)
    U[:,1] = V[:,1] / norm(V[:,1])
    for j in 2:N
        for jp in 1:(j-1)
            V[:,j] -= (J[:,j]'U[:,jp])*U[:,jp]
        end
        U[:,j] = V[:,j] / norm(V[:,j])
    end
    return U, V
end
A = rand(3,3)
B, C = gram_schmidt(A)
C
Q, R = qr(A)
# B[1,:]'B[2,:], B[2,:]'B[3,:], B[3,:]'B[1,:]
# norm(B[1,:]), norm(B[2,:]), norm(B[3,:])
# norm(B[:,1]), norm(B[:,2]), norm(B[:,3])

# C[:,1]'C[:,2], C[:,2]'C[:,3], C[:,3]'C[:,1]
# norm(C[1,:]), norm(C[2,:]), norm(C[3,:])
# norm(C[:,1]), norm(C[:,2]), norm(C[:,3])
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
    dt = 1e-3
    @showprogress @threads for dr = eachrow(schedules)[1:1:end]
        filename = "G:/DDM/lyapunov/lorenz/$(lpad(dr.idx, 5, '0')).csv"
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
    CSV.write("G:/DDM/lyapunov/lorenz.csv", result)
    plot(xticks = 0:20:100, legend = :none, size = [600, 300])
    plot!(result.ρ, result.λ1, lw = 2, color = 1)
    plot!(result.ρ, result.λ2, lw = 2, color = 2)
    plot!(result.ρ, result.λ3, lw = 2, color = 3)
    png("lyapunovGS")
end
lyapunov_lorenz()