include("ODEdata.jl")

const m = 10^(-3)
const μ = 10^(-6)
const R = 22 # 22
const L = 20m # 20m
const C = 47μ # 22μ
const T = 400μ
const γ = 11.7 # 11.75238
const η = 1309.5 # 1309.524
const RC = R*C

function factory_buck(idx::Int64, E::Number; flag_filesave = false)
    Vr(t) = γ + η * (mod(t, T))

    EdL = E/L
    
    controlterm = 0.0
    function buck(v::AbstractVector)
        V, I = v

        V̇ = - V/(RC) + I/C
        İ = - (V/L) + controlterm
        return [V̇, İ]
    end
    dt = 10^(-7); tend = 0.5
    t_ = 0:dt:tend
    Vr_ = Vr.(t_)

    ndatapoints = round(Int64, tend/(100dt))

    len_t_ = length(t_)
    traj = zeros(4, len_t_+1)
    u = [12.0, 0.55]
    du = buck(u)
    traj[1:2, 1] = u

    
    for t in 1:length(t_)
        controlterm = ifelse(u[1] < Vr_[t], EdL, 0)
        u, du = Euler(buck, u, dt)
        if t ≥ ndatapoints
            traj[3:4,   t] = du
            traj[1:2, t+1] = u
        end
    end
    traj = traj[:, 1:(end-1)]'

    if flag_filesave
        data = DataFrame(
            [t_ traj Vr_],
            ["t", "V", "I", "dV", "dI", "Vr"])
            @warn "file saving mode!"
            CSV.write("G:/buck/buck_$(lpad(idx, 6, '0')).csv", data)
        return nothing
    else
        data = DataFrame(
            [t_ traj Vr_][(end-ndatapoints):end, :],
            ["t", "V", "I", "dV", "dI", "Vr"])
    end

    return data
end

const κ = 400.0
const _μ = 172.363

function factory_soft(idx::Int64, d::Number)
    d2 = d/2
    
    function soft(tuv)
        t, u, v    = tuv

        impact = ifelse(abs(u) < d2, 0, -(κ^2)*sign(u)*(abs(u)-d2) - _μ*v)
        ṫ = 1
        u̇ = v
        v̇ = cospi(t) + impact
        return [ṫ, u̇, v̇]
    end
    dt = 10^(-5); tend = 50
    t_ = 0:dt:tend

    ndatapoints = round(Int64, tend/(10dt))

    len_t_ = length(t_)
    traj = zeros(6, len_t_+1)
    x = [0.0, 0.05853, 0.47898]
    dx = soft(x)
    traj[1:3, 1] = x

    
    for t in 1:length(t_)
        x, dx = RK4(soft, x, dt)
        if t ≥ ndatapoints
            traj[4:6,   t] = dx
            traj[1:3, t+1] = x
        end
    end
    traj = traj[:, (end-ndatapoints):(end-1)]'

    data = DataFrame(traj,
        ["t", "u", "v", "dt", "du", "dv"])

    return data
end