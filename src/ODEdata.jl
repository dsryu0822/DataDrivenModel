using CSV, DataFrames

function Euler(f::Function,v::AbstractVector, h=10^(-2))
    return v + h*f(v)
end
function RK4(f::Function,v::AbstractVector, h=10^(-2))
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4)
end


if !isfile("./data/buck.csv")
    m = 10^(-3)
    μ = 10^(-6)
    R = 22
    C = 47μ
    L = 20m
    T = 400μ
    γ = 11.75238
    η = 1309.524
    E = 53.500001

    Vr(t) = γ + η * (mod(t, T))
    global controlterm = 0
    function buck(v)
        V, I, _ = v

        V̇ = - V/(R*C) + I/C
        # İ = - (V/L) + ifelse(V < Vrt, E/L, 0)
        İ = - (V/L) + controlterm
        return [V̇, İ, 1]
    end

    u0 = [12.3, 0.55, 0.0]
    u_ = [u0]
    ∇_ = [buck(u_[end])]
    dt = 0.00001; tend = 0.25
    for t in dt:dt:tend
        # global Vrt = Vr(t)
        global controlterm = ifelse(u_[end][1] < Vr(t), E/L, 0)
        push!(u_, RK4(buck, u_[end], dt))
        push!(∇_, buck(u_[end]))
    end
    # push!(∇_, buck(u_[end]))
    U = stack(u_)[Not(end), Not(end)]
    ∇ = stack(∇_)[Not(end), Not(end)]
    CSV.write("data/buck.csv", DataFrame(
        [collect(dt:dt:tend)'; U; ∇; Vr.(dt:dt:tend)']'
      , ["t", "V", "I", "dV", "dI", "Vr"]))
end

if !isfile("./data/lorenz.csv")
    using OrdinaryDiffEq
    function lorenz(vec, p, t)
        x, y, z = vec

        ẋ = 10.0 * (y - x)
        ẏ = x * (28.0 - z) - y
        ż = x * y - (8 / 3) * z

        return [ẋ, ẏ, ż]
    end
    p,t = Nothing, Nothing
    u0 = [1; 0; 0]
    u0 = [-8; 8; 27]
    tspan = (0.0, 200.0)
    dt = 0.001
    prob = ODEProblem(lorenz, u0, tspan)
    sol = solve(prob, DP5(), saveat = dt, reltol = 1e-12, abstol = 1e-12)
    CSV.write("data/lorenz.csv", DataFrame(Array(sol)', :auto))
end