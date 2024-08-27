include("DDM.jl")

const _k = 1e+3
const _m = 1e-3
const _μ = 1e-6
const _n = 1e-9

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

# function factory_lorenz(ρ::Number; ic = [10.,10.,10.], tspan = [0., 10.], dt = 1e-4)
#     σ = 10
#     β = 8/3
#     function sys(v::AbstractVector)
#         x, y, z = v
#         dx = σ*(y - x)
#         dy = x*(ρ - z) - y
#         dz = x*y - β*z
#         return [dx, dy, dz]
#     end
        
#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)
    
#     t, tk = .0, 0
#     v = ic; DIM = length(v)
#     traj = zeros(len_t_+2, 2DIM)
#     while tk ≤ len_t_
#         x,y,z = v
#         v, dv = RK4(sys, v, dt)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  v
#             traj[tk  , DIM .+ (1:DIM)] = dv
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_lorenz(T::Type, args...; kargs...) =
# DataFrame(factory_lorenz(args...; kargs...), ["x", "y", "z", "dx", "dy", "dz"])

# const _R = 22 # 22
# const _L = 20_m # 20m
# const _C = 47_μ # 22μ
# const _T = 400μ
# const _γ = 11.7 # 11.75238
# const _η = 1309.5 # 1309.524
# const _RC = _R*_C
# Vr(t) = _γ + _η * (mod(t, _T))
# function factory_buck(E::Number; ic = [12.0, 0.55], tspan = [0.00, 0.01], dt = 1e-7)    
#     function sys(VI::AbstractVector, nonsmooth::Real)
#         V, I = VI

#         V̇ = - V/(_RC) + I/_C
#         İ = - (V/_L) + nonsmooth
#         return [V̇, İ]
#     end
#     EdL = E/_L
    
#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)

#     t, tk = .0, 0, 0
#     v = ic; DIM = length(v)
#     traj = zeros(len_t_+2, 2DIM)
#     while tk ≤ len_t_
#         V, I = v
#         nonsmooth = ifelse(V < Vr(t), EdL, 0); t += dt
#         v, dv = RK4(sys, v, dt, nonsmooth)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  v
#             traj[tk  , DIM .+ (1:DIM)] = dv
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_buck(T::Type, args...; kargs...) = 
# DataFrame(factory_buck(args...; kargs...), ["V", "I", "dV", "dI"])

function factory_soft(d::Number; ic = [.0, .05853, .47898], tspan = [0, 10], dt = 1e-5)
    κ = 400.0
    μ = 172.363
    function soft(tuv::AbstractVector, nonsmooth::Real)
        t, u, v = tuv
    
        ṫ = 1
        u̇ = v
        v̇ = cospi(t) + nonsmooth
        return [ṫ, u̇, v̇]
    end
    d2 = d/2

    t_ = first(tspan):dt:last(tspan)
    len_t_ = length(t_)
    
    t, tk = .0, 0
    x = ic; DIM = length(x)
    traj = zeros(len_t_+2, 2DIM)
    while tk ≤ len_t_
        t,u,v = x
        nonsmooth = ifelse(abs(u) < d2, 0, -(κ^2)*sign(u)*(abs(u)-d2) - μ*v)
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

function factory_hrnm(_f::Number; ic = [0.0, 0.0, 0.0, 0.1], tspan = [0, 100], dt = 1e-3)
    _a = 1.0
    _b = 3.0
    _c = 1.0
    _d = 5.0
    __k = 0.9
    # _f = 0.1
    _ω = 1.0
    _α = 0.1
    _β = 0.8 
    function sys(txyz::AbstractVector, nonsmooth::Real)
        t,x,y,z=txyz
    
        ṫ = 1
        ẋ = y - _a*x^3 + _b*x^2 + __k*x*z + _f*cos(_ω*t)
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
        v, dv = RK4(sys, v, dt, nonsmooth)

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

# function factory_gear(Fe::Number; ic = [0.1, 0.1, 0.1, 0.0], tspan = [0, 1000], dt = 1e-2)
#     __b = 1
#     ζ = 0.06
#     k1 = 0.06
#     ωh = 1
#     H = ωh^2
#     Fn = 0.3
#     φ = (√5 - 1)/2
#     function sys(xvΩθ::AbstractVector, nonsmooth::Real)
#         x,v,Ω,θ=xvΩθ
        
#         ẋ = v
#         v̇ = Fn + Fe*H*(cos(θ) + cos(Ω)) - (1 + k1*cos(Ω))*nonsmooth - 2*ζ*v
#         Ω̇ = ωh
#         θ̇ = φ
#         return [ẋ, v̇, Ω̇ , θ̇ ]
#     end
    
#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)
    
#     t, tk = .0, 0
#     u = ic; DIM = length(u)
#     traj = zeros(len_t_+2, 2DIM)
#     while tk ≤ len_t_
#         x,v,Ω,θ = u
#         t += dt
#         nonsmooth = ifelse(x > __b, x - __b, ifelse(x < -__b, x + __b, 0))
#         u, du = RK4(sys, u, dt, nonsmooth)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  u
#             traj[tk  , DIM .+ (1:DIM)] = du
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_gear(T::Type, args...; kargs...) = 
# DataFrame(factory_gear(args...; kargs...), ["x", "v", "Ω", "θ", "dx", "dv", "dΩ", "dθ"])

vrbl = [:dx1, :dx2, :dx3, :dx4], [:x1, :x2, :x3, :x4]
dt = 1e-2; tend = 10000;
function J_(x1, x2, x3, x4, β)
    L = 19.1_m  ; Ga1 = -0.0009302325
    C = 12.5_n  ; Ga2 = -0.000240577
    ν = 8300    ; # F = 0.3535533654213462

    F = 4.566348639778012
    G = √(C/(L*β))
    R = 1/G
    a₁ = Ga1*R
    a₂ = Ga2*R
    f = F*β
    ω = 2π*ν*C/G
    return [ 0                           1  0              0
             0 ifelse(abs(x1) ≤ 1, a₁, a₂)  1              0
             0                          -β -β -f*ω*cos(ω*x4)
             0                           0  0              0 ]
end

β = 0.163613
data = factory_mlcc(DataFrame, β, ic = [4.5, 4.5, 4.5, 0]; tspan = [0, tend], dt = 1e-3);
plot(data[900000:100:end,:].x1, data[900000:100:end,:].x2, color = :black, xlims = [-20, 140], xticks = -20:40:140)
vline!([-1, 1], color = :blue)
plot(data[900000:100:end,:].x1, data[900000:100:end,:].x2, color = :black, xlims = [-2, 3])
vline!([-1, 1], color = :blue)
@time λ = lyapunov_exponent(data[:, last(vrbl)], J_, β, T = tend)
# data = factory_mlcc(DataFrame, 0.339267; tspan = [0, 10000], dt = 0.013387)

data = factory_mlcc(DataFrame, .339267, ic = [0.0, 0.0, 0.0, 0]; tspan = [0, tend], dt);
plot(data[100000:100:end,:].x2, data[100000:100:end,:].x3)

function factory_mlcc(β::Number; ic = [1.0, -0.1, 0.1, 0.0], tspan = [0, 1000], dt = 1e-2)
    L = 19.1_m  ; Ga1 = -0.0009302325
    C = 12.5_n  ; Ga2 = -0.000240577
    ν = 8300    ; # F = 0.3535533654213462
    
    F = 4.566348639778012
    G = √(C/(L*β))
    R = 1/G
    a₁ = Ga1*R
    a₂ = Ga2*R
    f = F*β
    ω = 2π*ν*C/G
    # β = C/(L*(G^2))
    function sys(𝐱::AbstractVector, nonsmooth::Real)
        x₁,x₂,x₃,x₄=𝐱

        ẋ₁ = x₂
        ẋ₂ = x₃ - nonsmooth*x₂
        ẋ₃ = -β*(x₂ + x₃) + f*sin(ω*x₄)
        ẋ₄ = 1
        return [ẋ₁, ẋ₂, ẋ₃, ẋ₄]
    end
    
    t_ = first(tspan):dt:last(tspan)
    len_t_ = length(t_)
    
    t, tk = .0, 0
    v = ic; DIM = length(v)
    traj = zeros(len_t_+2, 2DIM)
    println("β = $β")
    # println("F = $F")
    # println("G = $G")
    # println("R = $R")
    println("a₁ = $a₁")
    println("a₂ = $a₂")
    println("f = $f")
    println("ω = $ω")
    # println("ν = $ν")
    while tk ≤ len_t_
        x₁,x₂,x₃,x₄ = v
        t += dt
        nonsmooth = ifelse(abs(x₁) ≤ 1, a₁, a₂)
        v, dv = RK4(sys, v, dt, nonsmooth)

        if t ≥ first(t_)
            tk += 1
            traj[tk+1,         1:DIM ] =  v
            traj[tk  , DIM .+ (1:DIM)] = dv
        end
    end
    return traj[2:(end-2), :]
end
factory_mlcc(T::Type, args...; kargs...) = 
DataFrame(factory_mlcc(args...; kargs...), ["x1", "x2", "x3", "x4", "dx1", "dx2", "dx3", "dx4"])