using CSV, DataFrames, ProgressMeter
@info "Packages CSV, DataFrames, ProgressMeter loaded"

# function Euler(f::Function, v::AbstractVector, h=10^(-2))
#     V1 = f(v)
#     return v + h*V1, V1
# end
function RK4(f, v::AbstractVector, h=10^(-2))
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end
function RK4(f::Function, v::AbstractVector, h=10^(-2), nonsmooth=0.0)
    V1 = f(v, nonsmooth)
    V2 = f(v + (h/2)*V1, nonsmooth)
    V3 = f(v + (h/2)*V2, nonsmooth)
    V4 = f(v + h*V3, nonsmooth)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end
function solve(f_, v, h = 10^(-2), t_ = nothing, DT = nothing, anc_ = nothing)
    bit_anc = anc_ |> isnothing
    V = zeros(1+length(t_), length(v))
    V[1, :] = v
    @showprogress for k in 1:length(t_)
        _v = bit_anc ? v : [v; anc_[k]]
        s = predict(DT, _v)
        v, _ = RK4(f_[s], v, h)
        V[k+1, :] = v
    end
    return V[1:(end-1), :]
end


const _m = 10^(-3)
const _μ = 10^(-6)
const _R = 22 # 22
const _L = 20_m # 20m
const _C = 47_μ # 22μ
const _T = 400_μ
const _γ = 11.7 # 11.75238
const _η = 1309.5 # 1309.524
const _RC = _R*_C
function buck(VI::AbstractVector, nonsmooth::Real)
    V, I = VI

    V̇ = - V/(_RC) + I/_C
    İ = - (V/_L) + nonsmooth
    return [V̇, İ]
end
Vr(t) = _γ + _η * (mod(t, _T))
function factory_buck(idx::Int64, E::Number, ic = [12.0, 0.55], tspan = (0.00, 0.01))
    EdL = E/_L
    
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
        x, dx = RK4(buck, x, dt, nonsmooth)
        if tk+1 ≥ (len_t_ - ndatapoints)
            traj[        1:DIM , tk+1] =  x
            traj[DIM .+ (1:DIM), tk  ] = dx
        end
    end
    traj = traj[:, 1:(end-1)]'
    traj = [t_ traj Vr_][(end-ndatapoints):end, :]

    return traj
end
factory_buck(T::Type, args...) = DataFrame(factory_buck(args...), ["t", "V", "I", "dV", "dI", "Vr"])

const __κ = 400.0
const __μ = 172.363
function soft(tuv::AbstractVector, nonsmooth::Real)
    t, u, v = tuv

    ṫ = 1
    u̇ = v
    v̇ = cospi(t) + nonsmooth
    return [ṫ, u̇, v̇]
end
function factory_soft(idx::Int64, d::Number, ic = [0.0, 0.0, 0.0], tspan = (0, 10))
    d2 = d/2
    
    dt = 10^(-5)
    t_ = 0:dt:last(tspan)
    ndatapoints = count(first(tspan) .< t_ .≤ last(tspan))
    len_t_ = length(t_)

    x = ic; DIM = length(x)
    traj = zeros(2DIM, len_t_+1)
    traj[1:DIM, 1] = x
    
    for tk in 1:len_t_
        nonsmooth = ifelse(abs(x[2]) < d2, 0, -(__κ^2)*sign(x[2])*(abs(x[2])-d2) - __μ*x[3])
        x, dx = RK4(soft, x, dt, nonsmooth)
        if tk+1 ≥ (len_t_ - ndatapoints)
            traj[        1:DIM , tk+1] =  x
            traj[DIM .+ (1:DIM), tk  ] = dx
        end
    end
    traj = traj[:, (end-ndatapoints-1):(end-2)]'

    return traj
end
factory_soft(T::Type, args...) = DataFrame(factory_soft(args...), ["t", "u", "v", "dt", "du", "dv"])

const _a = 1.0
const _b = 3.0
const _c = 1.0
const _d = 5.0
const _k = 0.9
const _f = 0.1
const _ω = 1.0
const _α = 0.1
const _β = 0.8 
function HR(txyz::AbstractVector, nonsmooth::Real)
    t,x,y,z=txyz

    ṫ = 1
    ẋ = y - _a*x^3 + _b*x^2 + _k*x*z + _f*cos(_ω*t)
    ẏ = _c - _d*x^2 - y
    ż = _α*nonsmooth + _β*x
    return [ṫ, ẋ, ẏ, ż]
end
function factory_HR(idx::Int64, a::Number, ic = [0.0, 0.1, 0.1, 0.1], tspan = (0, 20))
    dt = 10^(-3)
    t_ = 0:dt:last(tspan)

    ndatapoints = count(first(tspan) .< t_ .≤ last(tspan))

    len_t_ = length(t_)
    x = ic; DIM = length(x)
    traj = zeros(2DIM, len_t_+1)
    traj[1:DIM, 1] = x

    
    for tk in 1:len_t_
        nonsmooth = sign(x[4]+1) + sign(x[4]-1) - x[4]
        # if x[4] < -1
        #     nonsmooth = -2 -x[4]
        # elseif x[4] > 1
        #     nonsmooth = 2 - x[4]
        # end
        x, dx = RK4(HR, x, dt, nonsmooth)
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