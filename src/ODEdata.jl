using CSV, DataFrames, ProgressBars

# function Euler(f::Function, v::AbstractVector, h=10^(-2))
#     V1 = f(v)
#     return v + h*V1, V1
# end
function RK4(f::Function, v::AbstractVector, h=10^(-2))
    V1 = f(v)
    V2 = f(v + (h/2)*V1)
    V3 = f(v + (h/2)*V2)
    V4 = f(v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end
function RK4(f::Function, v::AbstractVector, nonsmooth::Real, h=10^(-2))
    V1 = f(v, nonsmooth)
    V2 = f(v + (h/2)*V1, nonsmooth)
    V3 = f(v + (h/2)*V2, nonsmooth)
    V4 = f(v + h*V3, nonsmooth)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end
function RK4(f::Function, s::Integer, v::AbstractVector, h=10^(-2))
    V1 = f(s, v)
    V2 = f(s, v + (h/2)*V1)
    V3 = f(s, v + (h/2)*V2)
    V4 = f(s, v + h*V3)
    return v + (h/6)*(V1 + 2V2 + 2V3 + V4), V1
end

# if !isfile("./data/lorenz.csv")
#     using OrdinaryDiffEq
#     function lorenz(vec, p, t)
#         x, y, z = vec

#         ẋ = 10.0 * (y - x)
#         ẏ = x * (28.0 - z) - y
#         ż = x * y - (8 / 3) * z

#         return [ẋ, ẏ, ż]
#     end
#     p,t = Nothing, Nothing
#     u0 = [1; 0; 0]
#     u0 = [-8; 8; 27]
#     tspan = (0.0, 200.0)
#     dt = 10^(-3)
#     prob = ODEProblem(lorenz, u0, tspan)
#     sol = solve(prob, DP5(), saveat = dt, reltol = 1e-12, abstol = 1e-12)
#     CSV.write("data/lorenz.csv", DataFrame(Array(sol)', :auto))
# end