using OrdinaryDiffEq
using CSV, DataFrames

function f(vec, p, t)
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
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, DP5(), saveat = dt, reltol = 1e-12, abstol = 1e-12)
CSV.write("test/lorenz.csv", DataFrame(Array(sol)', :auto))