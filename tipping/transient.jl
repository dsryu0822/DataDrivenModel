include("../core/header.jl")

"""''''''''''''''''''''''''''''''''''''''''''''''''''''

                    bifurcation diagram

''''''''''''''''''''''''''''''''''''''''''''''''''"""''
hrzn = []
vrtl = []
@async @showprogress for k in 0.9:0.0001:1
    for _ in 1:500
        traj = factory_foodchain(DataFrame, k, ic = [0, .85, .2 + .8rand(), .8], saveat = 0:1e-1:10000)[90000:end,:]
        bits = arglmin(traj.P)
        if maximum(traj.P) > 0.55
            push!(hrzn, fill(k, length(bits)))
            push!(vrtl, traj.P[bits])
            break
        else
            # println("k: $k should be tried again")
        end            
    end
end
bfdf = DataFrame(hrzn = [hrzn...;], vrtl = [vrtl...;])
# CSV.write("data_bifurcation_actual_foodchain.csv", DataFrame(hrzn = [hrzn...;], vrtl = [vrtl...;]))
bfdf = CSV.read("data_bifurcation_actual_foodchain.csv", DataFrame)
scatter(bfdf.hrzn, bfdf.vrtl, ylims = [0.55, 0.8], yticks = [0.55, 0.8], xticks = [0.9, .95, .96, 1], xlims = [0.9, 1], 
    ms = .5, msw = 0, color = :black, alpha = 0.5, legend = :none, ygrid = false,
    size = [400, 200])
png("bifurcaiton.png")


# """''''''''''''''''''''''''''''''''''''''''''''''''''''

#     transient lifetime distribution after critical point

# ''''''''''''''''''''''''''''''''''''''''''''''''''"""''
# tau_ = zeros(100000)
# @showprogress @threads for k in eachindex(tau_)
#     traj = factory_foodchain(DataFrame, 0.99976 + 2e-4, ic = rand(4), saveat = 0:1e-0:5000)
#     tau = findlast(diff(cumsum(traj.P .< 0.55)) .== 0)
#     tau_[k] = !isnothing(tau) ? tau : -1
# end
# freq = [count(x .≤ tau_ .< x + 250) for x in 250:250:4750]
# scatter(250:250:4750, freq ./ sum(freq), yscale = :log10, xticks = 0:1000:5000, xlims = [0, 5000], yticks = [1e-1, 1e-2], color = :white, legend = :none, xlabel = "transient lifetime", ylabel = "frequency")
# tau_[tau_ .== nothing] .= 5000
# tau_ = tau_[(tau_ .!= -1) .& (tau_ .< 5000)]
# histogram(tau_[100 .< tau_ .< 5000], yscale = :log10, bins = 25, normalize = true)

"""''''''''''''''''''''''''''''''''''''''''''''''''''''

                    transient chaos

''''''''''''''''''''''''''''''''''''''''''''''''''"""''
ic_chaos = [0.0, 0.9409608959542848, 0.38622912152951383, 0.7186582935906595]
traj1 = factory_foodchain(DataFrame, 0.99976 - 4e-5, ic = ic_chaos, saveat = 0:1e-1:10000)[1:10:end,:]
traj2 = factory_foodchain(DataFrame, 0.99976 + 4e-5, ic = ic_chaos, saveat = 0:1e-1:10000)[1:10:end,:]
teargs = (; color = :black, legend = :none, ylims = [0, 1.1], yticks = [0, 0.55, 1.1], xlims = [0, 9000])

plt_te_1 = plot(traj1.P; teargs...)
plt_te_2 = plot(traj2.P; teargs...)
plt_te = plot(plt_te_1, plt_te_2, layout = (2, 1), right_margin = 2mm, size = [400, 500])


"""''''''''''''''''''''''''''''''''''''''''''''''''''''

                    Recover by SINDy

''''''''''''''''''''''''''''''''''''''''''''''''''''"""
largs = (; xlabel = L"R", ylabel = L"C", zlabel = L"P", xticks = [0.3, 0.8], yticks = [0.2, 0.45], zticks = [0.65, 1.0])

traj0 = factory_foodchain(DataFrame, 0.95, ic = [0, 0.820915, 0.158239, 0.953786], saveat = 0:1e-1:5000)
traj1 = factory_foodchain(DataFrame, 0.96, ic = [0, 0.820915, 0.158239, 0.953786], saveat = 0:1e-1:5000)
_traj0 = traj0[traj0.t .≥ 4000,:]
_traj1 = traj1[traj1.t .≥ 4000,:]
plt_pp_1 = plot(_traj0.R, _traj0.C, _traj0.P; largs..., color = :black)
plt_pp_2 = plot(_traj1.R, _traj1.C, _traj1.P; largs..., color = :black)

vrbl = reverse(half(names(traj0)[2:end]))
cnfg = cook(vrbl; poly = 0:4)

@time f0 = SINDy(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print
@time f1 = SINDy(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print

tspan = (0, 5000)
prob0 = ODEProblem(define(Function, f0), collect(traj0[1, vrbl[2]]), tspan)
sol0 = solve(prob0);
prob1 = ODEProblem(define(Function, f1), collect(traj1[1, vrbl[2]]), tspan)
sol1 = solve(prob1);
plt_pp_3 = plot(eachrow(stack(sol0.u[1:20]))...; ticks = true, largs..., color = :red)
plt_pp_4 = plot(eachrow(stack(sol1.u[1:20]))...; ticks = true, largs..., color = :red)
# plot(
#     plot(eachrow(stack(sol0.u[1:10]))..., color = :red),
#     plot(eachrow(stack(sol1.u[1:10]))..., color = :red),
#     legend = :none, layout = (1, :)
# )



"""''''''''''''''''''''''''''''''''''''''''''''''''''

                Recover by implicit-SINDy

''''''''''''''''''''''''''''''''''''''''''''''''''"""
vrbl = reverse(half(names(traj0)[2:end]))
cnfg = cookPI(vrbl; poly = 0:3)

@time f0 = SINDyPI(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print
@time f1 = SINDyPI(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print

define(f0) |> print
define(f1) |> print

tspan = (0, 5000)
sargs = (; reltol = 1e-6, initializealg = DiffEqBase.BrownFullBasicInit())
prob0 = DAEProblem(define(Function, f0), collect(traj0[1, vrbl[1]]), collect(traj0[1, vrbl[2]]), tspan, differential_vars = ones(Bool, length(vrbl[1])))
sol0 = solve(prob0, Sundials.IDA(); sargs...);
prob1 = DAEProblem(define(Function, f1), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), tspan, differential_vars = ones(Bool, length(vrbl[1])))
sol1 = solve(prob1, Sundials.IDA(); sargs...);
plt_pp_5 = plot(eachrow(stack(sol0.u))...; largs..., color = :blue)
plt_pp_6 = plot(eachrow(stack(sol1.u))...; largs..., color = :blue)
# plot(
#     plot(eachrow(stack(sol0.u))...; largs..., color = :blue),
#     plot(eachrow(stack(sol1.u))...; largs..., color = :blue),
#     legend = :none, layout = (1, :)
# )
plot(plt_pp_1, plt_pp_2, plt_pp_3, plt_pp_4, plt_pp_5, plt_pp_6,
    layout = (3, 2), size = [600, 900], legend = :none)

"""''''''''''''''''''''''''''''''''''''''''''''''''''

            Bifucation Recover by implicit-SINDy

''''''''''''''''''''''''''''''''''''''''''''''''''"""
sargs = (; reltol = 1e-6, initializealg = DiffEqBase.BrownFullBasicInit(), saveat = 4000:1:5000)

# α = ((0.997) - 0.95)/0.01
α_ = -5:1e-2:5
hrzn = [[] for k in eachindex(α_)]
vrtl = [[] for k in eachindex(α_)]
f_ = []
for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg[1], method = "SINDyPI")

    push!(f_, define(Function, g; fname = "f_$k", sigdigits = 12))
end
@showprogress @threads for k in eachindex(α_)
    for l in 1:1000
        probg = DAEProblem(f_[k], rand(3), rand(3), (0, 5000), differential_vars = ones(Bool, length(vrbl[1])))
        solg = solve(probg, Sundials.IDA(); sargs...)
        trajg = DataFrame(stack(solg.u)', [:R, :C, :P])

        bits = arglmin(trajg.P)
        if maximum(trajg.P) > 0.55
            hrzn[k] = fill(k, length(bits))
            vrtl[k] = trajg.P[bits]
            break
        end
    end
end
bfdf = DataFrame(hrzn = [hrzn...;], vrtl = [vrtl...;])
# CSV.write("data_bifurcation_predict_foodchain.csv", bfdf)
bfdf = CSV.read("data_bifurcation_predict_foodchain.csv", DataFrame)
scatter(bfdf.hrzn, bfdf.vrtl, ms = .5, msw = 0, color = :blue, alpha = 0.5, legend = :none
, ylims = [0.55, 0.8], xformatter = x -> (x-500)/100, xticks = [0, 500, 600, 1000], yticks = [.55, .8], xlims = [0, 1000]
, size = [400, 200])
png("bifurcation.png")



"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Critical point

''''''''''''''''''''''''''''''''''''''''''''''''''"""
(
4.7778656471520655
+
4.777865647524595
) / 2
begin
α = 4.7778656471520655
probg = DAEProblem(define(Function, syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg[1], method = "SINDyPI"), sigdigits = 12), ones(3), [.85, .22, .8], (0, 10000), differential_vars = ones(Bool, length(vrbl[1])))
solg = solve(probg, Sundials.IDA(); sargs...)
trajg = DataFrame(stack(solg.u)', [:R, :C, :P])
plot(trajg.P, label = :none, color = :black, ylims = [0, 1.1], xlims = [0, 10000], xlabel = L"t", ylabel = L"P", title = L"α = %$α", margin = 5mm)
png("G:/temp.png")
end


"""''''''''''''''''''''''''''''''''''''''''''''''''''

        scaling law of transient lifetime for K

''''''''''''''''''''''''''''''''''''''''''''''''''"""
function food_chain(du, u, p, t)
    R, C, P = u
    K = p[1]
    xc, yc, xp, yp, R0, C0 = 0.4, 2.009, 0.08, 2.876, 0.16129, 0.5

    du[1] = R*(1 - frac(R,K)) - xc*yc*frac(C*R, R + R0)
    du[2] = xc*C*(frac(yc*R, R + R0) - 1) - xp*yp*frac(P*C, C + C0)
    du[3] = xp*P*(frac(yp*C, C + C0) - 1)
end
traj = factory_foodchain(DataFrame, .99976, ic = [0, .5, .5, .8], saveat = 0:1e-1:10000)[40000:end,:]
icc = shuffle([[r...] for r in eachrow(traj[:, [:R, :C, :P]])])
ε_ = exp10.(range(-5, -2, 20))
tau2__ = [zeros(1000) for k in eachindex(ε_)]
for k in eachindex(ε_)
    @showprogress @threads for l in eachindex(tau2__[k])
        prob = ODEProblem(food_chain, icc[l], (0, 100000), [0.99976+ε_[k]])
        sol = solve(prob; sargs...)
        tau = findfirst(sol[3, :] .< 0.55)
        tau2__[k][l] = !isnothing(tau) ? sol.t[tau] : last(sol.t)
    end
end
plt_scl2 = scatter(ε_, mean.(tau2__), scale = :log10, size = [400, 400], color = :white, legend = :none, xticks = exp10.(-10:-1))



"""''''''''''''''''''''''''''''''''''''''''''''''''''

        scaling law of transient lifetime for α

''''''''''''''''''''''''''''''''''''''''''''''''''"""
α_c = 4.777865647
# probg = DAEProblem(define(Function, syntheticSINDy((1-α_c)*f0.matrix + α_c*f1.matrix, vrbl, cnfg[1], method = "SINDyPI"), sigdigits = 12),
#                       ones(3), [.85, .21, .8], (0, 10000), differential_vars = ones(Bool, length(vrbl[1])))
# solg = solve(probg, Sundials.IDA(); sargs...)
# trajg = DataFrame(stack(solg.u)', [:R, :C, :P])
# icc = shuffle([[r...] for r in eachrow(trajg)])
# define(String, syntheticSINDy((1-α_c)*f0.matrix + α_c*f1.matrix, vrbl, cnfg[1], method = "SINDyPI"), sigdigits = 12)

ε_ = exp10.(range(-5, -2, 20))
α_ = α_c .+ ε_
f_ = []
for k in eachindex(α_)
    α = α_[k]
    g = syntheticSINDy((1-α)*f0.matrix + α*f1.matrix, vrbl, cnfg[1], method = "SINDyPI")
    push!(f_, define(Function, g; fname = "f_$k", sigdigits = 16))
end
tau1__ = [zeros(10) for k in eachindex(α_)]
for k in eachindex(α_)
    @showprogress @threads for l in eachindex(tau1__[k])
        probg = DAEProblem(f_[k], ones(3), icc[l], (0, 1e+7), differential_vars = ones(Bool, length(vrbl[1])))
        solg = solve(probg, Sundials.IDA(); sargs...)
        t_ = solg.t
        u_ = solg.u
        for _ in 1:100
            if last(u_)[3] > .55
                probg = DAEProblem(f_[k], ones(3), icc[l], (last(t_), 1e+7), differential_vars = ones(Bool, length(vrbl[1])))
                solg = solve(probg, Sundials.IDA(); sargs...)
                t_ = solg.t
                u_ = solg.u
            end
        end
        tau = findfirst(stack(u_)[3, :] .< 0.55)
        tau1__[k][l] = !isnothing(tau) ? t_[tau] : last(t_)
    end
end
# plt_scl1 = scatter(ε_, mean.(tau1__), scale = :log10, size = [400, 400], color = :white, msc = :blue, legend = :none, xticks = exp10.(-10:-1))

# JLD2.@save "temp2.jld2" ε_ tau1__ tau2__
# JLD2.@load "temp1.jld2"
# maximum.(tau1__)
plt_trs = plot(size = [400, 400], xticks = exp10.(-10:-1), ylims = [1e+1, 1e+6]);
scatter!(plt_trs, ε_, mean.(tau2__), scale = :log10, color = :white, msc = :black, shape = :circ,  label = "actual");
scatter!(plt_trs, ε_, mean.(tau1__), scale = :log10, color = :white, msc = :blue,  shape = :rect, label = "recovered");
plt_trs


tau2__[1] |> histogram
begin
    include("../core/header.jl")
    include("../core/factorio.jl")
    using DifferentialEquations
    import Sundials
    import DiffEqBase

    traj0 = factory_foodchain(DataFrame, 0.95, ic = [0, 0.820915, 0.158239, 0.953786], saveat = 0:1e-1:5000)
    traj1 = factory_foodchain(DataFrame, 0.96, ic = [0, 0.820915, 0.158239, 0.953786], saveat = 0:1e-1:5000)

    vrbl = reverse(half(names(traj0)[2:end]))
    cnfg = cookPI(vrbl; poly = 0:3)

    @time f0 = SINDyPI(traj0, vrbl, cnfg, λ = 1e-8); f0 |> print
    @time f1 = SINDyPI(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print

    condition(u,t,integrator) = 0.1 - u[3]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition,affect!)
    sargs = (; reltol = 1e-24, initializealg = DiffEqBase.BrownFullBasicInit(), callback = cb)
end
