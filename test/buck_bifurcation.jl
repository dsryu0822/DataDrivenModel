const m = 10^(-3)
const μ = 10^(-6)
const R = 22 # 22
const L = 20m # 20m
const C = 47μ # 22μ
const T = 400μ
const γ = 11.7 # 11.75238
const η = 1309.5 # 1309.524
const RC = R*C

using ProgressBars
packages = [:DataFrames, :CSV, :Plots, :LaTeXStrings]
@time for package in ProgressBar(packages)
    @eval using $(package)
end
println(join(packages, ", "), " loaded!")

function Euler(f::Function,v::AbstractVector, h=10^(-2))
    V1 = f(v)
    return v + h*V1, V1
end
function RK4(f::Function,v::AbstractVector, h=10^(-2))
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end

Vr(t) = γ + η * (mod(t, T))

function factory_buck(idx::Int64, E::Number)
    EdL = E/L
    
    controlterm = 0.0
    function buck(v::AbstractVector)
        V, I = v

        V̇ = - V/(RC) + I/C
        İ = - (V/L) + controlterm
        # İ = - (V/L) + ifelse(V < Vr(t), E/L, 0)
        return [V̇, İ]
    end
    dt = 10^(-7); tend = 0.5
    t_ = 0:dt:tend
    Vr_ = Vr.(t_)
    ndatapoints = round(Int64, 1/200dt)

    len_t_ = length(t_)
    traj = zeros(4, len_t_+1)
    u = [12.0, 0.55]
    du = buck(u)
    traj[1, 1:2] = u

    
    for t in 1:length(t_)
        controlterm = ifelse(u[1] < Vr_[t], EdL, 0)
        u, du = Euler(buck, u, dt)
        if t ≥ ndatapoints
            traj[3:4,   t] = du
            traj[1:2, t+1] = u
        end
    end
    traj = traj[:, 1:(end-1)]'
    
    data = DataFrame(
        [t_ traj Vr_][(end-ndatapoints):end, :],
        ["t", "V", "I", "dV", "dI", "Vr"])
    # CSV.write("G:/buck/buck_$(lpad(idx, 6, '0')).csv", data)
    return data
end

# E_range = 15:40
E_range = [15:0.1:24; 24:0.01:40]
schedule = DataFrame(idx = eachindex(E_range), E = E_range)

xdots = Float64[]; ydots = Float64[]
for dr in ProgressBar(eachrow(schedule))
    data = factory_buck(dr.idx, dr.E)
    
    # idx_jump = (diff(data.dI)) .< -100
    # sampledV = data.V[Not([end])][idx_jump]
    idx_sampled = diff(data.Vr) .< 0
    sampledV = data[Not(1), :V][idx_sampled]
    append!(xdots, fill(dr.E, length(sampledV)))
    append!(ydots, sampledV)
end
@showtime a1 = scatter(xdots, ydots, label = :none, msw = 0, color = :black, ms = 1, alpha = 0.5, size = (1000, 600))

png(a1, "a1")
plot(a1)

# schedule[100:800,:]
# dr = eachrow(schedule)[92]
# # @time data = factory_buck(dr.idx, dr.E); # to check compilation time
# @time data = factory_buck(dr.idx, dr.E)
# plot(data.V, data.I, alpha = 0.9, legend = :best, label = :none, xlabel = L"V", ylabel = L"I", color = :black, size = (400, 400), title = L"E = 30.5")
# png("a2")

# idx_sampled = diff(data.Vr) .< 0
# sampledV = data[Not(1), :V][idx_sampled]
# scatter!(data[Not(1), :V][idx_sampled], data[Not(1), :I][idx_sampled], alpha = 0.9, label = "Sampled", color = :black, ms = 1, msw = 0)
# # idx_jump = (diff(data.dI)) .< -100;
# # scatter!(data.V[Not(end)][idx_jump], data.I[Not(end)][idx_jump], alpha = 0.9)