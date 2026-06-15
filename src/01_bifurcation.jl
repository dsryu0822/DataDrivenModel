include("../core/header.jl")

sol = factory_lorenz63(DataFrame, [10, 28, 8/3])[200000:end, :]
plot(sol.x, sol.y, sol.z, alpha = .5)

σ_ = range(6, 15, length = 1001)
ρ_ = range(120, 150, length = 1001)
β_ = range(3, 5, length = 1001)
if !isfile("bifurcation_lorenz63.jld2")
     @info "Calculating bifurcation data for Lorenz63..."
    bfcn = Dict{Float64, Vector{Float64}}()
    @showprogress @threads for k in eachindex(σ_)
        sol = factory_lorenz63(DataFrame, [σ_[k], ρ_[k], β_[k]], ic = [100, 100, 100])
        z_ = sol.z[sol.t .≥ 900]
        bfcn[k] = z_[arglmax(z_)]
    end
    JLD2.@save "bifurcation_lorenz63.jld2" bfcn
else
    @info "Loading bifurcation data for Lorenz63 from file..."
    JLD2.@load "bifurcation_lorenz63.jld2" bfcn
end

scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .1, msw = 0, color = :black); png("temp")

traj0 = factory_lorenz63(DataFrame, [σ_[1], ρ_[1], β_[1]], tspan = 900:1e-2:1000)
traj1 = factory_lorenz63(DataFrame, [σ_[101], ρ_[101], β_[101]], tspan = 900:1e-2:1000)

vrbl = reverse(half(names(traj0[:, Not(:t)])))
cnfg = cook(last(vrbl), poly = 0:2)
f0 = SINDy(traj0, vrbl, cnfg; λ = 1e-3); f0 |> print
f1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); f1 |> print

rcvd0 = ssolve(f0, rand(3), 900:1e-3:1000)
rcvd1 = ssolve(f1, rand(3), 900:1e-3:1000)
plot(
    plot(traj0.x, traj0.y, traj0.z, color = :black),
    plot(traj1.x, traj1.y, traj1.z, color = :black),
    plot(rcvd0[:, 2], rcvd0[:, 3], rcvd0[:, 4], color = :blue),
    plot(rcvd1[:, 2], rcvd1[:, 3], rcvd1[:, 4], color = :blue),
)

α_ = 0:1e-2:10
f_ = [syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg, method = "SINDy") for α in α_]
bfcn = Dict{Float64, Vector{Float64}}()
@showprogress @threads for k in eachindex(α_)
    try
        sol = solve(ODEProblem(define(Function, f_[k]; fname = "f_$k", sigdigits = 20), [100, 100, 100], (0, 1000)), saveat = 900:1e-3:1000)
        z_ = sol[3, :]
        bfcn[α_[k]] = z_[arglmax(z_)]
    catch e
        continue
    end
end
JLD2.@save "bifurcation_lorenz63_recover.jld2" bfcn
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .1, msw = 0, color = :blue)

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                        lorenz96

''''''''''''''''''''''''''''''''''''''''''''''''''"""
function factory_lorenz96(A::Number; ic = rand(6), tspan = 0:1e-3:2500)
    function lorenz96(du, u, p, t)
        A = p[1]

        for k in eachindex(ic)
            du[k] = u[mod1(k-1, 6)]*(u[mod1(k+1, 6)] - u[mod1(k-2, 6)]) - u[k] + A
        end
        return du
    end
    # sol = solve(ODEProblem(lorenz96, ic, (0, last(tspan)), [A]), RK4(), dt = tspan.step.hi, adaptive=false)
    sol = solve(ODEProblem(lorenz96, ic, (0, last(tspan)), [A]), saveat = tspan)
    # sol = solve(ODEProblem(lorenz96, ic, (0, last(tspan)), [A]), reltol = 1e-24)
    matrix = Matrix([sol.t'; sol[:, :]; stack([lorenz96(zeros(6), u, [A], 0) for u in sol.u])]')
    return matrix[sol.t .≥ first(tspan), :][1:end-1, :]
end
factory_lorenz96(T::Type, args...; kargs...) =
DataFrame(factory_lorenz96(args...; kargs...), ["t", "x1", "x2", "x3", "x4", "x5", "x6", "dx1", "dx2", "dx3", "dx4", "dx5", "dx6"])

sol = factory_lorenz96(DataFrame, 6)[200000:100:end, :]
plot(sol.x1, sol.x2, sol.x3, alpha = .1)

if !isfile("bifurcation_lorenz96.jld2")
     @info "Calculating bifurcation data for Lorenz96..."
     bfcn = Dict{Float64, Vector{Float64}}()
     @showprogress @threads for A in 5.3:1e-3:6
         sol = factory_lorenz96(DataFrame, A)
         x_ = sol.x1[sol.t .≥ 2000]
         bfcn[A] = x_[arglmax(x_)]
     end
     JLD2.@save "bifurcation_lorenz96.jld2" bfcn
 else
     @info "Loading bifurcation data for Lorenz96 from file..."
     JLD2.@load "bifurcation_lorenz96.jld2" bfcn
 end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .1, msw = 0, color = :black, xticks = 5.3:0.1:6)
png("bifurcation_lorenz96.png")

traj0 = factory_lorenz96(DataFrame, 5.3, tspan = 2000:1e-2:2500)
traj1 = factory_lorenz96(DataFrame, 5.4, tspan = 2000:1e-2:2500)

vrbl = reverse(half(names(traj0[:, Not(:t)])))
cnfg = cook(last(vrbl), poly = 0:2)
f0 = SINDy(traj0, vrbl, cnfg; λ = 1e-3); f0 |> print
f1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); f1 |> print

rcvd0 = ssolve(f0, rand(6), 2000:1e-3:2500)
rcvd1 = ssolve(f1, rand(6), 2000:1e-3:2500)
plot(
    plot(traj0.x1, traj0.x2, traj0.x3, color = :black),
    plot(traj1.x1, traj1.x2, traj1.x3, color = :black),
    plot(rcvd0[:, 2], rcvd0[:, 3], rcvd0[:, 4], color = :blue),
    plot(rcvd1[:, 2], rcvd1[:, 3], rcvd1[:, 4], color = :blue),
)
plot(traj0.x1, xlims = [0, 1000])
plot(traj1.x1, xlims = [0, 1000])

α_ = 0:1e-2:7
f_ = [syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg, method = "SINDy") for α in α_]
bfcn = Dict{Float64, Vector{Float64}}()
@showprogress @threads for k in eachindex(α_)
    try
        sol = solve(ODEProblem(define(Function, f_[k]; fname = "f_$k", sigdigits = 20), rand(6), (0, 2500)), saveat = 2000:1e-3:2500)
        x_ = sol[2, :]
        bfcn[α_[k]] = x_[arglmax(x_)]
    catch e
        continue
    end
end
JLD2.@save "bifurcation_lorenz96_recover.jld2" bfcn

filter(x -> !isempty(second(x)), bfcn)
for (k, v) in bfcn
    if v |> isempty
        delete!(bfcn, k)
    end
end
scatter(dict2bifurcation(bfcn)..., ms = .5, ma = .1, msw = 0, color = :blue, xticks = 0:7); png("bifurcation_lorenz96_recover.png")

asdf = rand(10)
plot(asdf, ticks = 4, lims = 4)