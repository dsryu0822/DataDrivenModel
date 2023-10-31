include("ODEdata.jl")

const _m = 10^(-3)
const _μ = 10^(-6)
const _R = 22 # 22
const _L = 20_m # 20m
const _C = 47_μ # 22μ
const _T = 400_μ
const _γ = 11.7 # 11.75238
const _η = 1309.5 # 1309.524
const _RC = _R*_C

function factory_buck(idx::Int64, E::Number, tspan = (0.495, 0.5), ic = [12.0, 0.55])
    Vr(t) = _γ + _η * (mod(t, _T))

    EdL = E/_L
    
    function buck(VI::AbstractVector, nonsmooth::Real)
        V, I = VI

        V̇ = - V/(_RC) + I/_C
        İ = - (V/_L) + nonsmooth
        return [V̇, İ]
    end
    dt = 10^(-7)
    t_ = 0:dt:last(tspan)
    Vr_ = Vr.(t_)

    ndatapoints = count(first(tspan) .< t_ .≤ last(tspan))

    len_t_ = length(t_)
    x = ic; DIM = length(x)
    traj = zeros(2DIM, len_t_+1)
    traj[1:DIM, 1] = x
    
    for tk in 1:length(t_)
        nonsmooth = ifelse(x[1] < Vr_[tk], EdL, 0)
        x, dx = RK4(buck, x, nonsmooth, dt)
        if tk+1 ≥ (len_t_ - ndatapoints)
            traj[        1:DIM , tk+1] =  x
            traj[DIM .+ (1:DIM), tk  ] = dx
        end
    end
    traj = traj[:, 1:(end-1)]'
    traj = [t_ traj Vr_][(end-ndatapoints):end, :]

    data = DataFrame(traj,
        ["t", "V", "I", "dV", "dI", "Vr"])

    return data
end

const __κ = 400.0
const __μ = 172.363
function factory_soft(idx::Int64, d::Number, tspan = (45, 50), ic = [0.0, 0.05853, 0.47898])
    d2 = d/2
    
    function soft(tuv::AbstractVector, nonsmooth::Real)
        t, u, v = tuv

        ṫ = 1
        u̇ = v
        v̇ = cospi(t) + nonsmooth
        return [ṫ, u̇, v̇]
    end
    dt = 10^(-5)
    t_ = 0:dt:last(tspan)

    ndatapoints = count(first(tspan) .< t_ .≤ last(tspan))

    len_t_ = length(t_)
    x = ic; DIM = length(x)
    traj = zeros(2DIM, len_t_+1)
    traj[1:DIM, 1] = x

    
    for tk in 1:len_t_
         # u == x[2], v == x[3]
        nonsmooth = ifelse(abs(x[2]) < d2, 0, -(__κ^2)*sign(x[2])*(abs(x[2])-d2) - __μ*x[3])
        x, dx = RK4(soft, x, nonsmooth, dt)
        if tk+1 ≥ (len_t_ - ndatapoints)
            traj[        1:DIM , tk+1] =  x
            traj[DIM .+ (1:DIM), tk  ] = dx
        end
    end
    traj = traj[:, (end-ndatapoints):(end-1)]'

    data = DataFrame(traj,
        ["t", "u", "v", "dt", "du", "dv"])

    return data
end

const _a = 1.0
const _b = 3.0
const _c = 1.0
const _d = 5.0
const _k = 0.9
const _f = 0.1
const _ω = 1.0
const _α = 0.1
const _β = 0.8
function factory_HR(idx::Int64, a::Number, tspan = (90, 100), ic = [0.0, 0.1, 0.1, 0.1])
    
    function HR(txyz::AbstractVector, nonsmooth::Real)
        t,x,y,z=txyz

        ṫ = 1
        ẋ = y - _a*x^3 + _b*x^2 + _k*x*z + _f*cos(_ω*t)
        ẏ = _c - _d*x^2 - y
        ż = _α*nonsmooth + _β*x
        return [ṫ, ẋ, ẏ, ż]
    end
    dt = 10^(-3)
    t_ = 0:dt:last(tspan)

    ndatapoints = count(first(tspan) .< t_ .≤ last(tspan))

    len_t_ = length(t_)
    x = ic; DIM = length(x)
    traj = zeros(2DIM, len_t_+1)
    traj[1:DIM, 1] = x

    
    for tk in 1:len_t_
         # u == x[2], v == x[3]
        nonsmooth = sign(x[4]+1) + sign(x[4]-1) - x[4]
        # if x[4] < -1
        #     nonsmooth = -2 -x[4]
        # elseif x[4] > 1
        #     nonsmooth = 2 - x[4]
        # end
        x, dx = RK4(HR, x, nonsmooth, dt)
        if tk+1 ≥ (len_t_ - ndatapoints)
            traj[        1:DIM , tk+1] =  x
            traj[DIM .+ (1:DIM), tk  ] = dx
        end
    end
    traj = traj[:, (end-ndatapoints):(end-1)]'

    data = DataFrame(traj,
        ["t", "x", "y", "z", "dt", "dx", "dy", "dz"])

    return data
end