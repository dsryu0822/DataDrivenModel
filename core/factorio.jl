include("DDM.jl")
function euler(f::Function, v::AbstractVector, h=1e-2)
    V1 = f(v)
    return v + h*V1, V1
end
function RK4(J::AbstractMatrix, U::AbstractMatrix, dt=1e-2)
    V1 = J*U
    V2 = J*(U + (dt/2)*V1)
    V3 = J*(U + (dt/2)*V2)
    V4 = J*(U + dt*V3)
    return U + (dt/6)*(V1 + 2V2 + 2V3 + V4)
end
function RK4(f::Union{Function, STLSQresult}, v::AbstractVector, h=1e-2, nonsmooth=nothing)
    if nonsmooth |> isnothing
        V1 = f(v)
        V2 = f(v + (h/2)*V1)
        V3 = f(v + (h/2)*V2)
        V4 = f(v + h*V3)
    else
        V1 = f(v, nonsmooth)
        V2 = f(v + (h/2)*V1, nonsmooth)
        V3 = f(v + (h/2)*V2, nonsmooth)
        V4 = f(v + h*V3, nonsmooth)
    end
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end

function solve(f_, v, h = 1e-2, t_ = nothing, DT = nothing, anc_ = nothing)
    V = zeros(length(t_), length(v))
    for k in eachindex(t_)
        V[k, :] = v
        s = apply_tree(DT, [v; anc_[k]])
        v, _ = RK4(f_[s], v, h)
    end
    return V
end
function solve(f_, v, h = 1e-2, t_ = nothing, DT = nothing)
    V = zeros(length(t_), length(v))
    for k in eachindex(t_)
        V[k, :] = v
        s = apply_tree(DT, v)
        v, _ = RK4(f_[s], v, h)
    end
    return V
end


function factory_lorenz(ρ::Number; ic = [10.,10.,10.], tspan = [0., 10.], dt = 1e-4)
    σ = 10
    β = 8/3
    function lorenz(v::AbstractVector)
        x, y, z = v
        dx = σ*(y - x)
        dy = x*(ρ - z) - y
        dz = x*y - β*z
        return [dx, dy, dz]
    end

    t_ = first(tspan):dt:last(tspan)
    if first(tspan) |> !iszero
        ic[1] = first(tspan)
    end
    
    ndatapoints = count(first(tspan) .≤ t_ .< last(tspan))
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
factory_lorenz(T::Type, args...; kargs...) =
DataFrame(factory_lorenz(args...; kargs...), ["x", "y", "z", "dx", "dy", "dz"])

const _m = 1e-3
const _μ = 1e-6
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
function factory_buck(E::Number; ic = [12.0, 0.55], tspan = [0.00, 0.01], dt = 1e-7)
    EdL = E/_L
    
    t_ = first(tspan):dt:last(tspan)
    len_t_ = length(t_)

    t, tk = .0, 0, 0
    v = ic; DIM = length(v)
    traj = zeros(len_t_+2, 2DIM)
    while tk ≤ len_t_
        V, I = v
        nonsmooth = ifelse(V < Vr(t), EdL, 0); t += dt
        v, dv = RK4(buck, v, dt, nonsmooth)

        if t ≥ first(t_)
            tk += 1
            traj[tk+1,         1:DIM ] =  v
            traj[tk  , DIM .+ (1:DIM)] = dv
        end
    end
    return traj[2:(end-2), :]
end
factory_buck(T::Type, args...; kargs...) = 
DataFrame(factory_buck(args...; kargs...), ["V", "I", "dV", "dI"])

const __κ = 400.0
const __μ = 172.363
function soft(tuv::AbstractVector, nonsmooth::Real)
    t, u, v = tuv

    ṫ = 1
    u̇ = v
    v̇ = cospi(t) + nonsmooth
    return [ṫ, u̇, v̇]
end
function factory_soft(d::Number; ic = [.0, .05853, .47898], tspan = [0, 10], dt = 1e-5)
    d2 = d/2

    t_ = first(tspan):dt:last(tspan)
    len_t_ = length(t_)
    
    t, tk = .0, 0
    x = ic; DIM = length(x)
    traj = zeros(len_t_+2, 2DIM)
    while tk ≤ len_t_
        t,u,v = x
        nonsmooth = ifelse(abs(u) < d2, 0, -(__κ^2)*sign(u)*(abs(u)-d2) - __μ*v)
        x, dx = RK4(soft, x, dt, nonsmooth)

        if t ≥ first(t_)
            tk += 1
            traj[tk+1,         1:DIM ] =  x
            traj[tk  , DIM .+ (1:DIM)] = dx
        end
    end
    return traj[2:(end-2), :]
end
factory_soft(T::Type, args...; kargs...) =
DataFrame(factory_soft(args...; kargs...), ["t", "u", "v", "dt", "du", "dv"])

const _a = 1.0
const _b = 3.0
const _c = 1.0
const _d = 5.0
const _k = 0.9
# const _f = 0.1
const _ω = 1.0
const _α = 0.1
const _β = 0.8 
function factory_hrnm(_f::Number; ic = [0.0, 0.0, 0.0, 0.1], tspan = [0, 100], dt = 1e-3)
    function HR(txyz::AbstractVector, nonsmooth::Real)
        t,x,y,z=txyz
    
        ṫ = 1
        ẋ = y - _a*x^3 + _b*x^2 + _k*x*z + _f*cos(_ω*t)
        ẏ = _c - _d*x^2 - y
        ż = _α*nonsmooth + _β*x
        return [ṫ, ẋ, ẏ, ż]
    end
    
    t_ = first(tspan):dt:last(tspan)
    len_t_ = length(t_)
    
    t, tk = .0, 0
    v = ic; DIM = length(v)
    traj = zeros(len_t_+2, 2DIM)
    while tk ≤ len_t_
        t,x,y,z = v
        nonsmooth = sign(z+1) + sign(z-1) - z
        v, dv = RK4(HR, v, dt, nonsmooth)

        if t ≥ first(t_)
            tk += 1
            traj[tk+1,         1:DIM ] =  v
            traj[tk  , DIM .+ (1:DIM)] = dv
        end
    end
    return traj[2:(end-2), :]
end
factory_hrnm(T::Type, args...; kargs...) = 
DataFrame(factory_hrnm(args...; kargs...), ["t", "x", "y", "z", "dt", "dx", "dy", "dz"])
