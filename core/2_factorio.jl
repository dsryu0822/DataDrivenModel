# include("DDM.jl")

# const _k = 1e+3
# const _m = 1e-3
# const _μ = 1e-6
# const _n = 1e-9

# function euler(f::Function, v::AbstractVector, h=1e-2)
#     V1 = f(v)
#     return v + h*V1, V1
# end
# function RK4(J::AbstractMatrix, U::AbstractMatrix, h=1e-2)
#     V1 = J*U
#     V2 = J*(U + (h/2)*V1)
#     V3 = J*(U + (h/2)*V2)
#     V4 = J*(U + h*V3)
#     return U + (h/6)*(V1 + 2V2 + 2V3 + V4)
# end
# function RK4(f::Function, v::AbstractVector, h=1e-2)
#     V1 = f(v)
#     V2 = f(v + (h/2)*V1)
#     V3 = f(v + (h/2)*V2)
#     V4 = f(v + h*V3)
#     return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
# end
# function RK4(f::Union{Function, STLSQresult}, v::AbstractVector, h=1e-2, nonsmooth=nothing)
#     if nonsmooth |> isnothing
#         V1 = f(v)
#         V2 = f(v + (h/2)*V1)
#         V3 = f(v + (h/2)*V2)
#         V4 = f(v + h*V3)
#     else
#         V1 = f(v, nonsmooth)
#         V2 = f(v + (h/2)*V1, nonsmooth)
#         V3 = f(v + (h/2)*V2, nonsmooth)
#         V4 = f(v + h*V3, nonsmooth)
#     end
#     return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
# end

# function solve(f_, v, h = 1e-2, t_ = nothing, DT = nothing, anc_ = nothing)
#     V = zeros(length(t_), length(v))
#     for k in eachindex(t_)
#         V[k, :] = v
#         s = apply_tree(DT, [v; anc_[k]])
#         v, _ = RK4(f_[s], v, h)
#     end
#     return V
# end
# function solve(f_, v, t_, DT = nothing)
#     apply_ = typeof(DT) <: Root ? apply_tree : apply_forest
#     h = t_[2] - t_[1]
#     V = zeros(length(t_), length(v))
#     if DT |> isnothing
#         for k in eachindex(t_)
#             V[k, :] = v
#             v, _ = RK4(f_, v, h)
#         end
#     else
#         for k in eachindex(t_)
#             V[k, :] = v
#             s = apply_(DT, v)
#             v, _ = RK4(f_[s], v, h)
#         end
#     end
#     return V
# end

# function factory_soft(d::Number; ic = [-1, .05853, -.47898], tspan = [0, 10], dt = 1e-5)
#     κ = 400.0
#     μ = 172.363
#     function soft(tuv::AbstractVector, nonsmooth::Real)
#         t, u, v = tuv
    
#         ṫ = 1
#         u̇ = v
#         v̇ = cospi(t) + nonsmooth
#         return [ṫ, u̇, v̇]
#     end
#     d2 = d/2

#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)
    
#     t, tk = .0, 0
#     x = ic; DIM = length(x)
#     traj = zeros(len_t_+2, 2DIM)
#     while tk ≤ len_t_
#         t,u,v = x
#         nonsmooth = ifelse(abs(u) < d2, 0, -(κ^2)*sign(u)*(abs(u)-d2) - μ*v)
#         x, dx = RK4(soft, x, dt, nonsmooth)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  x
#             traj[tk  , DIM .+ (1:DIM)] = dx
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_soft(T::Type, args...; kargs...) =
# DataFrame(factory_soft(args...; kargs...), ["t", "u", "v", "dt", "du", "dv"])

# function factory_hrnm(_f::Number; ic = [-1, 0.0, 0.0, 0.1], tspan = [0, 100], dt = 1e-3)
#     _a = 1.0
#     _b = 3.0
#     _c = 1.0
#     _d = 5.0
#     _k = 0.9
#     # _f = 0.1
#     _ω = 1.0
#     _α = 0.1
#     _β = 0.8 
#     function sys(txyz::AbstractVector, nonsmooth::Real)
#         t,x,y,z=txyz
    
#         ṫ = 1
#         ẋ = y - _a*x^3 + _b*x^2 + _k*x*z + _f*cos(_ω*t)
#         ẏ = _c - _d*x^2 - y
#         ż = _α*nonsmooth + _β*x
#         return [ṫ, ẋ, ẏ, ż]
#     end
    
#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)
    
#     t, tk = .0, 0
#     v = ic; DIM = length(v)
#     traj = zeros(len_t_+2, 2DIM)
#     while tk ≤ len_t_
#         t,x,y,z = v
#         nonsmooth = sign(z+1) + sign(z-1) - z
#         v, dv = RK4(sys, v, dt, nonsmooth)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  v
#             traj[tk  , DIM .+ (1:DIM)] = dv
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_hrnm(T::Type, args...; kargs...) = 
# DataFrame(factory_hrnm(args...; kargs...), ["t", "x", "y", "z", "dt", "dx", "dy", "dz"])

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

##########################################################################
#                                                                        #
#                               not used                                 #
#                                                                        #
##########################################################################

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

# function factory_mlcc(β::Number; ic = [1.0, -0.1, 0.1, 0.0], tspan = [0, 1000], dt = 1e-2)
#     L = 19.1_m  ; Ga1 = -0.0009302325
#     C = 12.5_n  ; Ga2 = -0.000240577
#     ν = 8300    ; # F = 0.3535533654213462
    
#     F = 4.566348639778012
#     G = √(C/(L*β))
#     R = 1/G
#     a₁ = Ga1*R
#     a₂ = Ga2*R
#     f = F*β
#     ω = 2π*ν*C/G
#     # β = C/(L*(G^2))
#     function sys(𝐱::AbstractVector, nonsmooth::Real)
#         x₁,x₂,x₃,x₄=𝐱

#         ẋ₁ = x₂
#         ẋ₂ = x₃ - nonsmooth*x₂
#         ẋ₃ = -β*(x₂ + x₃) + f*sin(ω*x₄)
#         ẋ₄ = 1
#         return [ẋ₁, ẋ₂, ẋ₃, ẋ₄]
#     end
    
#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)
    
#     t, tk = .0, 0
#     v = ic; DIM = length(v)
#     traj = zeros(len_t_+2, 2DIM)

#     while tk ≤ len_t_
#         x₁,x₂,x₃,x₄ = v
#         t += dt
#         nonsmooth = ifelse(abs(x₁) ≤ 1, a₁, a₂)
#         v, dv = RK4(sys, v, dt, nonsmooth)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  v
#             traj[tk  , DIM .+ (1:DIM)] = dv
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_mlcc(T::Type, args...; kargs...) = 
# DataFrame(factory_mlcc(args...; kargs...), ["x1", "x2", "x3", "x4", "dx1", "dx2", "dx3", "dx4"])

# function factory_epid(ε::Number; ic = [.0, .0331, .0001, .9668], tspan = [0, 100], dt = 1e-4)
#     σ = 0.01
#     μ = 0.01
#     γ = 50
#     β₀ = 1510
#     T = 1
#     function sys(v::AbstractVector, nonsmooth::Real)
#         t,S,I,R=v

#         ṫ = 1
#         Ṡ = σ - μ*S - nonsmooth*S*I
#         İ = nonsmooth*S*I - (γ + μ)*I
#         Ṙ = γ*I - μ*R
#         return [ṫ, Ṡ, İ, Ṙ]
#     end
    
#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)
    
#     t, tk = .0, 0
#     v = ic; DIM = length(v)
#     traj = zeros(len_t_+2, 2DIM)
#     while tk ≤ len_t_
#         t,S,I,R=v
#         nonsmooth = β₀*(1 + ε*ifelse(mod(t, T) < .5, -1, 1))
#         v, dv = RK4(sys, v, dt, nonsmooth)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  v
#             traj[tk  , DIM .+ (1:DIM)] = dv
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_epid(T::Type, args...; kargs...) = 
# DataFrame(factory_epid(args...; kargs...), ["t", "S", "I", "R", "dt", "dS", "dI", "dR"])

# function factory_rkvt(a₁::Number; ic = [.5, .5], tspan = [0, 500], dt = 1e-3)
#     β = 0.5
#     d = 0.05
#     a₂ = 0.05
#     C₁ = 0.2
#     C₂ = 0.2
#     function sys(v::AbstractVector, nonsmooth::Real)
#         x₁,x₂=v

#         ẋ₁ = x₁*(1 - x₁ - x₂) - ifelse(x₁ > C₁, a₁, 0)
#         ẋ₂ = x₂*(β*x₁ - d) - ifelse(x₂ > C₂, a₂, 0)
#         return [ẋ₁, ẋ₂]
#     end
    
#     t_ = first(tspan):dt:last(tspan)
#     len_t_ = length(t_)
    
#     t, tk = .0, 0
#     v = ic; DIM = length(v)
#     traj = zeros(len_t_+2, 2DIM)
#     while tk ≤ len_t_
#         x₁,x₂=v
#         nonsmooth = 0
#         v, dv = RK4(sys, v, dt, nonsmooth)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  v
#             traj[tk  , DIM .+ (1:DIM)] = dv
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_rkvt(T::Type, args...; kargs...) = 
# DataFrame(factory_rkvt(args...; kargs...), ["x1", "x2", "dx1", "dx2"])


# function factory_buck(E::Number; ic = [12.0, 0.55], tspan = [0.00, 0.01], dt = 1e-7)    
#     μ = 1e-6
#     m = 1e-3
#     _R = 22 # 22
#     _L = 20_m # 20m
#     _C = 47_μ # 22μ
#     _T = 400μ
#     _γ = 11.7 # 11.75238
#     _η = 1309.5 # 1309.524
#     _RC = _R*_C
#     Vr(t) = _γ + _η * (mod(t, _T))
    
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
#     traj = zeros(len_t_+2, 2DIM+2)
#     while tk ≤ len_t_
#         V, I = v
#         Vrt = Vr(t)
#         nonsmooth = ifelse(V < Vrt, EdL, 0); t += dt
#         v, dv = RK4(sys, v, dt, nonsmooth)

#         if t ≥ first(t_)
#             tk += 1
#             traj[tk+1,         1:DIM ] =   v
#             traj[tk  , DIM .+ (1:DIM)] =  dv
#             traj[tk  ,       2DIM + 1] =   t
#             traj[tk  ,       2DIM + 2] = Vrt
#         end
#     end
#     return traj[2:(end-2), :]
# end
# factory_buck(T::Type, args...; kargs...) = 
# DataFrame(factory_buck(args...; kargs...), ["V", "I", "dV", "dI", "t", "Vr"])


# function factory_hastingpowell(a0::Number; ic = [0., 1, 1, 1], saveat = 0:1e-2:10)
#       b0, w0, d0, w1, d1,  w2, d2, a1, w3, d3,  c3 = (
#     0.05,  1, 10,  2, 10, 1.5, 10,  1,  2, 20, 0.7 )
#     function sys(v::AbstractVector)
#         _,v1,v2,v3 = v

#         dt = 1
#         dv1 = a0*v1 - b0*v1^2 - w0*frac(v1*v2, d0 + v1)
#         dv2 = w1*frac(v1*v2, d1 + v1) - w2*frac(v2*v3, d2 + v2) - a1*v2
#         dv3 = w3*frac(v2*v3, d3 + v2) - c3*v3
#         return [dt, dv1, dv2, dv3]
#     end
        
#     len_t_ = length(tspan); h = tspan.step.hi
#     v = ic; DIM = length(v)
#     traj = zeros(2+len_t_, 2DIM); traj[1, 1:DIM] = v

#     tk = 0
#     while tk ≤ len_t_
#         t,_,_,_ = v
#         v, dv = RK4(sys, v, h)

#         if t ≥ first(saveat)
#             tk += 1
#             traj[tk+1,         1:DIM ] =  v
#             traj[tk  , DIM .+ (1:DIM)] = dv
#         end
#     end
#     return traj[1:(end-2), :]
# end
# factory_hastingpowell(T::Type, args...; kargs...) =
# DataFrame(factory_hastingpowell(args...; kargs...), ["t", "v1", "v2", "v3", "dt", "dv1", "dv2", "dv3"])[:, Not(:dt)]

load_packages([:DataFrames, :DifferentialEquations, :OrdinaryDiffEqLowOrderRK, :DiffEqBase])

frac(x, y) = x ./ y

"""
    factory_lorenz63(σ, ρ, β; ic = [5, 5, 5], saveat = 0:1e-3:1000)

The Lorenz system is a system of ordinary differential equations first studied by Edward Lorenz.
It is notable for having chaotic solutions for certain parameter values and initial conditions.
"""
function factory_lorenz63(pp; ic = [5, 5, 5], saveat = 0:1e-3:1000)
    function sys(du, u, p, t)
        x, y, z = u; σ, ρ, β = p
        
        du[1] = σ*(y - x)
        du[2] = x*(ρ - z) - y
        du[3] = x*y - β*z
        return du
    end
    sol = solve(ODEProblem(sys, ic, (0, last(saveat)), pp); saveat)
    matrix = Matrix([sol.t'; sol[:, :]; stack([sys(zeros(3), u, pp, 0) for u in sol.u])]')
    return matrix[sol.t .≥ first(saveat), :][1:end-1, :]
end
factory_lorenz63(T::Type, args...; kargs...) =
DataFrame(factory_lorenz63(args...; kargs...), ["t", "x", "y", "z", "dx", "dy", "dz"])

function factory_foodchain(K::Number; ic = [.85, .12 + .34rand(), .8], saveat = 0:1e-2:10)
     xc,    yc,   xp,    yp,      R0,  C0 = (
    0.4, 2.009, 0.08, 2.876, 0.16129, 0.5)
    function sys(du, u, p, t)
        R,C,P = u
        K = p[1]

        du[1] = R*(1 - frac(R,K)) - xc*yc*frac(C*R, R + R0)
        du[2] = xc*C*(frac(yc*R, R + R0) - 1) - xp*yp*frac(P*C, C + C0)
        du[3] = xp*P*(frac(yp*C, C + C0) - 1)
        return du
    end
    sol = solve(ODEProblem(sys, ic, (0, last(saveat)), [K]), RK4(), dt = saveat.step.hi, adaptive=false, maxiters = 1e+7)
    matrix = Matrix([sol.t'; sol[:, :]; stack([sys(zeros(3), u, [K], 0) for u in sol.u])]')
    return matrix[sol.t .≥ first(saveat), :]
end
factory_foodchain(T::Type, args...; kargs...) =
DataFrame(factory_foodchain(args...; kargs...), ["t", "R", "C", "P", "dR", "dC", "dP"])

function factory_thomas(b::Number; ic = rand(3), saveat = 0:1e-2:2000)
    function sys(du, u, p, t)
        x, y, z = u
        b = p[1]

        du[1] = sin(y) - b*x
        du[2] = sin(z) - b*y
        du[3] = sin(x) - b*z
        return du
    end
    sol = solve(ODEProblem(sys, ic, (0, last(saveat)), [b]), RK4(), dt = saveat.step.hi, adaptive=false, maxiters = 1e+7)
    matrix = Matrix([sol.t'; sol[:, :]; stack([sys(zeros(3), u, [b], 0) for u in sol.u])]')
    return matrix[sol.t .≥ first(saveat), :][1:end-1, :]
end
factory_thomas(T::Type, args...; kargs...) =
DataFrame(factory_thomas(args...; kargs...), ["t", "x", "y", "z", "dx", "dy", "dz"])

function factory_pollinator(κ::Number; ic = [1.,1.], saveat = 0:1e-2:10)
      α, β,  _h,   γp,   γA = (
    0.3, 1, 0.4, 1.93, 1.77 )
    function sys(du, u, p, t)
        P,A = u
        κ = p[1]

        du[1] = P*(α - β*P + frac(γp*A, 1 + _h*γp*A))
        du[2] = A*(α - κ - β*A + frac(γA*P, 1 + _h*γA*P))
        return du
    end
    sol = solve(ODEProblem(sys, ic, (0, last(saveat)), [κ]), RK4(), dt = saveat.step.hi, adaptive=false, maxiters = 1e+7)
    matrix = Matrix([sol.t'; sol[:, :]; stack([sys(zeros(2), u, [κ], 0) for u in sol.u])]')
    return matrix[sol.t .≥ first(saveat), :][1:end-1, :]
end
factory_pollinator(T::Type, args...; kargs...) =
DataFrame(factory_pollinator(args...; kargs...), ["t", "P", "A", "dP", "dA"])

function factory_algaezooplankton(F::Number; ic = [5, 5], saveat = 0:1e-2:100)
      r, K,   g,  hA,   i,   e,    m,  hZ = (
    0.5, 6, 0.4, 0.6, 0.1, 0.6, 0.15, 0.5 )
    function sys(du, u, p, t)
        A,Z = u
        F = p[1]

        du[1] = r*A*(1 - frac(A,K)) - g*Z*frac(A, A + hA) + i*(K - A)
        du[2] = e*g*Z*frac(A, A + hA) - m*Z - F*frac(Z^2, Z^2 + hZ^2)
        return du
    end
    sol = solve(ODEProblem(sys, ic, (0, last(saveat)), [F]), RK4(), dt = saveat.step.hi, adaptive=false, maxiters = 1e+7)
    matrix = Matrix([sol.t'; sol[:, :]; stack([sys(zeros(2), u, [F], 0) for u in sol.u])]')
    return matrix[sol.t .≥ first(saveat), :][1:end-1, :]
end
factory_algaezooplankton(T::Type, args...; kargs...) =
DataFrame(factory_algaezooplankton(args...; kargs...), ["t", "A", "Z", "dA", "dZ"])
