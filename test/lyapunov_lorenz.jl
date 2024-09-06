using CSV, DataFrames, LinearAlgebra, ProgressMeter, Plots, Dates; @info now()
using Base.Threads: @threads # Base.Threads.nthreads()

function RK4(f::Function, v::AbstractVector, h=1e-3)
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end
function RK4(J::AbstractMatrix, U::AbstractMatrix, dt=1e-3)
    V1 = J*U
    V2 = J*(U + (dt/2)*V1)
    V3 = J*(U + (dt/2)*V2)
    V4 = J*(U + dt*V3)
    return U + (dt/6)*(V1 + 2V2 + 2V3 + V4)
end
function factory_lorenz(ρ::Number; ic = [10.,10.,10.], tspan = [0., 1000.], dt = 1e-3)
    σ = 10
    β = 8/3
    function sys(v::AbstractVector)
        x, y, z = v
        dx = σ*(y - x)
        dy = x*(ρ - z) - y
        dz = x*y - β*z
        return [dx, dy, dz]
    end
        
    t_ = first(tspan):dt:last(tspan)
    len_t_ = length(t_)
    
    t, tk = .0, 0
    v = ic; DIM = length(v)
    traj = zeros(len_t_+2, 2DIM)
    while tk ≤ len_t_
        v, dv = RK4(sys, v, dt)

        if t ≥ first(t_)
            tk += 1
            traj[tk+1,         1:DIM ] =  v
            traj[tk  , DIM .+ (1:DIM)] = dv
        end
    end
    return traj[2:(end-2), :]
end
factory_lorenz(T::Type, args...; kargs...) =
DataFrame(factory_lorenz(args...; kargs...), ["x", "y", "z", "dx", "dy", "dz"])

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
J_lorenz(x,y,z,ρ) = [
     -10   10  0
      ρ-z -1  -x
        y  x  -8/3
]

const dt = 1e-3
ρ_range = 0:.1:100
cd("G:/DDM")

function lyapunov_lorenz_GS()
    schedules = DataFrame(idx = eachindex(ρ_range), ρ = ρ_range)
    schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; 
    @showprogress @threads for dr = eachrow(schedules)
        data = factory_lorenz(DataFrame, dr.ρ; dt)
        
        λ = zeros(3);
        J = J_lorenz(collect(data[1, 1:3])..., dr.ρ)
        U = I(3)
        for i in 2:nrow(data)
            U, V = gram_schmidt(U)
            λ += V |> eachcol .|> norm .|> log          

            U = RK4(J, U, dt)
            J = J_lorenz(collect(data[i, 1:3])..., dr.ρ)
        end

        λ ./= dt*nrow(data)
        dr[[:λ1, :λ2, :λ3]] .= sort(λ, rev = true)
    end
    sort!(schedules, :ρ)
    CSV.write("lyapunov/lorenzGS.csv", schedules, bom = true)
    plot(xticks = 0:20:100, legend = :none, xlabel = "ρ", ylabel = "λ", size = [600, 300])
    plot!(schedules.ρ, schedules.λ1, lw = 2, color = 1)
    plot!(schedules.ρ, schedules.λ2, lw = 2, color = 2)
    plot!(schedules.ρ, schedules.λ3, lw = 2, color = 3)
    png("lyapunov/lorenzGS.png")
end

function lyapunov_lorenz_QR()
    schedules = DataFrame(idx = eachindex(ρ_range), ρ = ρ_range)
    schedules[!, :λ1] .= .0; schedules[!, :λ2] .= .0; schedules[!, :λ3] .= .0; 
    @showprogress @threads for dr = eachrow(schedules)
        data = factory_lorenz(DataFrame, dr.ρ; dt)
        
        λ = zeros(3);
        J = J_lorenz(collect(data[1, 1:3])..., dr.ρ)
        Q, _ = qr(J)
        for i in 2:nrow(data)
            Q, R = qr(Q)
            λ += R .|> abs |> diag .|> log   

            Q = RK4(J, Matrix(Q), dt)
            J = J_lorenz(collect(data[i, 1:3])..., dr.ρ)
        end

        λ ./= dt*nrow(data)
        dr[[:λ1, :λ2, :λ3]] .= sort(λ, rev = true)
    end
    sort!(schedules, :ρ)
    CSV.write("lyapunov/lorenzQR.csv", schedules, bom = true)
    plot(xticks = 0:20:100, legend = :none, xlabel = "ρ", ylabel = "λ", size = [600, 300])
    plot!(schedules.ρ, schedules.λ1, lw = 2, color = 1)
    plot!(schedules.ρ, schedules.λ2, lw = 2, color = 2)
    plot!(schedules.ρ, schedules.λ3, lw = 2, color = 3)
    png("lyapunov/lorenzQR.png")
end

lyapunov_lorenz_GS()
lyapunov_lorenz_QR()