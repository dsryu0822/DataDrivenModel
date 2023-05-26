using CSV, DataFrames
using Plots


Oₘ = CSV.read("G:/DDM/atopy.txt", DataFrame, delim = '\t', header = ["P","B","D","R","K","G","kappa","t"]) #, limit = 10000)
dt = Oₘ.t[2] - Oₘ.t[1]
Oₘ.B .*= 100

# Oₘ.P ./= sqrt(sum(abs2,Oₘ.P))
# Oₘ.B ./= sqrt(sum(abs2,Oₘ.B))
# Oₘ.D ./= sqrt(sum(abs2,Oₘ.D))


# plot(
#     plot(Oₘ.t[50000:60000], Oₘ.P[50000:60000], title = "P"),
#     plot(Oₘ.t[50000:60000], Oₘ.B[50000:60000], title = "B"),
#     plot(Oₘ.t[50000:60000], Oₘ.D[50000:60000], title = "D"),
# )

X = Matrix(Oₘ[1:1000:(end-1),1:3])
dXdt = (Matrix(Oₘ[2:1000:end,1:3]) - Matrix(Oₘ[1:1000:(end-1),1:3])) / dt
# plot(sum.(abs2, eachrow(diff(dXdt, dims = 1))))

Θ = libararize(X, 4)

using Plots

m = 10^(-3)
μ = 10^(-6)

R = 22
C = 4.7μ
L = 20m
T = 400μ
# γ = 11.75238
# η = 1309.524
γ = 11
η = 1300
E = 50

Vr(t) = γ + η*(mod(t, T))

function buck(v)
    V, I, t = v
    V̇ = -V/(R*C) + (I/C)
    İ = -(V/L) + ifelse(V ≥ Vr(t), 0, E/L)
    ṫ = 1
    return  [V̇,İ,ṫ]
end

dt = 0.000001
v0 = Float64[12.28, 0.55, 0]; v_ = [v0]
c_ = []
d_ = []
for t in 0:dt:0.0025
    push!(v_, RK4(buck, v_[end], dt))
    push!(c_, Vr(t))
    push!(d_, buck(v_[end]))
end
p1 = plot(stack(v_)[1,:], xformatter = (x -> dt*x), label = "V", title = "Time Evolution")
plot!(c_, label = "Vr", legend = :bottomright)
plot(stack(v_)[2,:], xformatter = (x -> dt*x))
plot(stack(d_)[2,:])

p2 = plot(stack(v_)[1,:], stack(v_)[2,:], title = "Phase Plane", xlabel = "V", ylabel = "I", legend = :none, color = :black)
plot(p1, p2, size = (800, 800), layout = @layout [a{0.3h}; b])

diff(stack(v_)[2,:])

extrema(stack(v_)[2,:])