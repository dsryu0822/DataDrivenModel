using OrdinaryDiffEq

function lorenz(u, p, t)
    x, y = u

    c = 0.5*exp((x-32000)/8000) #저항계수
    F = -c*y*abs(y) #저항력
    g = 9.8/(1+(-x+32000)/6370000)^2 #중력가속도

    ẋ = y
    ẏ = (F-70*g)/70 #운동방정식
    return [ẋ, ẏ]
end

u0 = [0, 32000]
tspan = (0.0, 200.0)
dt = 0.1
prob = ODEProblem(lorenz, u0, tspan)
# sol = solve(prob, Tsit5(), dense = true) # sol(sol.t,Val{4}); sol(sol.t,Val{5})
sol = solve(prob, RK4(), saveat = dt)

t_b=findfirst(ysol -> ysol<=0,sol[2,:])
sol[1:198]
print(t_b)