include.("../core/" .* readdir("core")[[1,2,3,4,6]])
cd("figure")

trajargs = (; ticks = [], color = :black)
rcvdargs = (; ticks = [], color = :blue)
failargs = (; ticks = [], color = :red)
bfcnargs = (; yticks = [], msw = 0, ma = .5, ms = .5)

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Lorenz-63

''''''''''''''''''''''''''''''''''''''''''''''''''"""
traj1 = factory_lorenz63(DataFrame, [6, 120, 3], saveat = 1900:1e-3:2000); CSV.write("traj1.csv", traj1);
traj2 = factory_lorenz63(DataFrame, [6.9, 123, 3.2], saveat = 1900:1e-3:2000); CSV.write("traj2.csv", traj2);
traj3 = factory_lorenz63(DataFrame, [12.3, 141, 4.4], saveat = 1900:1e-3:2000); CSV.write("traj3.csv", traj3);
vrbl = reverse(half(names(traj1[:, Not(:t)])))
cnfg = cook(vrbl, poly = 0:2)
f1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); f1 |> print
f2 = SINDy(traj2, vrbl, cnfg; λ = 1e-3); f2 |> print
f3 = SINDy(traj3, vrbl, cnfg; λ = 1e-3); f3 |> print
rcvd1 = solve(ODEProblem(define(Function, f1), collect(traj1[1, vrbl[2]]), (0, 1000)), maxiters = 1e+7, RK4(), dt = 1e-3, adaptive=false, saveat = 900:1e-3:1000) |> stack; CSV.write("rcvd1.csv", DataFrame(rcvd1', :auto));
rcvd2 = solve(ODEProblem(define(Function, f2), collect(traj2[1, vrbl[2]]), (0, 1000)), maxiters = 1e+7, RK4(), dt = 1e-3, adaptive=false, saveat = 900:1e-3:1000) |> stack; CSV.write("rcvd2.csv", DataFrame(rcvd2', :auto));
rcvd3 = solve(ODEProblem(define(Function, f3), collect(traj3[1, vrbl[2]]), (0, 1000)), maxiters = 1e+7, RK4(), dt = 1e-3, adaptive=false, saveat = 900:1e-3:1000) |> stack; CSV.write("rcvd3.csv", DataFrame(rcvd3', :auto));
plot(
    plot(eachcol(traj3[:, last(vrbl)])...; trajargs...),
    plot(eachrow(rcvd1)...; rcvdargs...),
    plot(eachrow(rcvd2)...; rcvdargs...),
    plot(eachrow(rcvd3)...; rcvdargs...)
)
png("trajectory_lorenz63.png")

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                        food chain

''''''''''''''''''''''''''''''''''''''''''''''''''"""
traj1 = factory_foodchain(DataFrame, 0.95, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000); CSV.write("traj1.csv", traj1);
traj2 = factory_foodchain(DataFrame, 0.96, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000); CSV.write("traj2.csv", traj2);
traj3 = factory_foodchain(DataFrame, 0.99, ic = [0.820915, 0.158239, 0.953786], saveat = 4000:1e-1:5000); CSV.write("traj3.csv", traj3);
vrbl = reverse(half(names(traj1)[2:end]))
cnfg = cookPI(vrbl; poly = 0:3)
f1 = SINDyPI(traj1, vrbl, cnfg, λ = 1e-8); f1 |> print
f2 = SINDyPI(traj2, vrbl, cnfg, λ = 1e-8); f2 |> print
f3 = SINDyPI(traj3, vrbl, cnfg, λ = 1e-8); f3 |> print
rcvd1 = solve(DAEProblem(define(Function, f1), collect(traj1[1, vrbl[1]]), collect(traj1[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1]))), Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6) |> stack; CSV.write("rcvd1.csv", DataFrame(rcvd1', :auto));
rcvd2 = solve(DAEProblem(define(Function, f2), collect(traj2[1, vrbl[1]]), collect(traj2[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1]))), Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6) |> stack; CSV.write("rcvd2.csv", DataFrame(rcvd2', :auto));
rcvd3 = solve(DAEProblem(define(Function, f3), collect(traj3[1, vrbl[1]]), collect(traj3[1, vrbl[2]]), (0, 5000), differential_vars = ones(Bool, length(vrbl[1]))), Sundials.IDA(), initializealg = DiffEqBase.BrownFullBasicInit(), reltol = 1e-6) |> stack; CSV.write("rcvd3.csv", DataFrame(rcvd3', :auto));
cnfg = cook(vrbl; poly = 0:4)
g1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); g1 |> print
g2 = SINDy(traj2, vrbl, cnfg; λ = 1e-3); g2 |> print
g3 = SINDy(traj3, vrbl, cnfg; λ = 1e-3); g3 |> print
fail1 = solve(ODEProblem(define(Function, g1), collect(traj1[1, vrbl[2]]), (0, 5000)), maxiters = 1e+7, RK4(), dt = 1e-1, adaptive=false, saveat = 4000:1e-1:5000) |> stack; CSV.write("fail1.csv", DataFrame(fail1', :auto));
fail2 = solve(ODEProblem(define(Function, g2), collect(traj2[1, vrbl[2]]), (0, 5000)), maxiters = 1e+7, RK4(), dt = 1e-1, adaptive=false, saveat = 4000:1e-1:5000) |> stack; CSV.write("fail2.csv", DataFrame(fail2', :auto));
fail3 = solve(ODEProblem(define(Function, g3), collect(traj3[1, vrbl[2]]), (0, 5000)), maxiters = 1e+7, RK4(), dt = 1e-1, adaptive=false, saveat = 4000:1e-1:5000) |> stack; CSV.write("fail3.csv", DataFrame(fail3', :auto));
plot(
    plot(eachcol(traj1[:, last(vrbl)])...; trajargs...),
    plot(eachcol(traj2[:, last(vrbl)])...; trajargs...),
    plot(eachcol(traj3[:, last(vrbl)])...; trajargs...),
    plot(eachrow(rcvd1)...; rcvdargs...),
    plot(eachrow(rcvd2)...; rcvdargs...),
    plot(eachrow(rcvd3)...; rcvdargs...),
    plot(eachrow(fail1)...; failargs...),
    plot(eachrow(fail2)...; failargs...),
    plot(eachrow(fail3)...; failargs...),
    size = (600, 600)
)
png("trajectory_foodchain.png")

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                        Thomas

''''''''''''''''''''''''''''''''''''''''''''''''''"""
ic0 = [0.0, 0.0, .1]
traj1 = factory_thomas(DataFrame, .140, ic = ic0, saveat = 1000:1e-2:2000); CSV.write("traj1.csv", traj1);
traj2 = factory_thomas(DataFrame, .141, ic = ic0, saveat = 1000:1e-2:2000); CSV.write("traj2.csv", traj2);
traj3 = factory_thomas(DataFrame, .150, ic = ic0, saveat = 1000:1e-2:2000); CSV.write("traj3.csv", traj3);
vrbl = reverse(half(names(traj1)[2:end]))
cnfg = cook(vrbl; poly = 0:1, trig = [1], format = sin)
f1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); f1 |> print
f2 = SINDy(traj2, vrbl, cnfg; λ = 1e-3); f2 |> print
f3 = SINDy(traj3, vrbl, cnfg; λ = 1e-3); f3 |> print
rcvd1 = solve(ODEProblem(define(Function, f1), collect(traj1[1, last(vrbl)]), (0, 2000)), maxiters = 1e+7, RK4(), dt = 1e-2, adaptive=false, saveat = 1000:1e-2:2000) |> stack; CSV.write("rcvd1.csv", DataFrame(rcvd1', :auto));
rcvd2 = solve(ODEProblem(define(Function, f2), collect(traj2[1, last(vrbl)]), (0, 2000)), maxiters = 1e+7, RK4(), dt = 1e-2, adaptive=false, saveat = 1000:1e-2:2000) |> stack; CSV.write("rcvd2.csv", DataFrame(rcvd2', :auto));
rcvd3 = solve(ODEProblem(define(Function, f3), collect(traj3[1, last(vrbl)]), (0, 2000)), maxiters = 1e+7, RK4(), dt = 1e-2, adaptive=false, saveat = 1000:1e-2:2000) |> stack; CSV.write("rcvd3.csv", DataFrame(rcvd3', :auto));
cnfg = cook(vrbl; poly = 0:5)
g1 = SINDy(traj1, vrbl, cnfg; λ = 1e-3); g1 |> print
g2 = SINDy(traj2, vrbl, cnfg; λ = 1e-3); g2 |> print
g3 = SINDy(traj3, vrbl, cnfg; λ = 1e-3); g3 |> print
fail1 = solve(ODEProblem(define(Function, g1), collect(traj1[1, vrbl[2]]), (0, 2000)), maxiters = 1e+7, RK4(), dt = 1e-2, adaptive=false, saveat = 1000:1e-2:2000) |> stack; CSV.write("fail1.csv", DataFrame(fail1', :auto));
fail2 = solve(ODEProblem(define(Function, g2), collect(traj2[1, vrbl[2]]), (0, 2000)), maxiters = 1e+7, RK4(), dt = 1e-2, adaptive=false, saveat = 1000:1e-2:2000) |> stack; CSV.write("fail2.csv", DataFrame(fail2', :auto));
fail3 = solve(ODEProblem(define(Function, g3), collect(traj3[1, vrbl[2]]), (0, 2000)), maxiters = 1e+7, RK4(), dt = 1e-2, adaptive=false, saveat = 1000:1e-2:2000) |> stack; CSV.write("fail3.csv", DataFrame(fail3', :auto));
plot(
    plot(eachcol(traj1[:, last(vrbl)])...; trajargs...),
    plot(eachcol(traj2[:, last(vrbl)])...; trajargs...),
    plot(eachcol(traj3[:, last(vrbl)])...; trajargs...),
    plot(eachrow(rcvd1)...; rcvdargs...),
    plot(eachrow(rcvd2)...; rcvdargs...),
    plot(eachrow(rcvd3)...; rcvdargs...),
    plot(eachrow(fail1)...; failargs...),
    plot(eachrow(fail2)...; failargs...),
    plot(eachrow(fail3)...; failargs...),
    size = (600, 600)
)
png("trajectory_thomas.png")


bfcn_lorenz63 = JLD2.load("bifurcation_lorenz63.jld2")["bfcn"]; CSV.write("bfcn.csv", DataFrame(stack(dict2bifurcation(bfcn_lorenz63)), ["h", "v"]))
bfcn_lorenz63_rcvd = JLD2.load("bifurcation_lorenz63_recover.jld2")["bfcn"]; CSV.write("bfcn_.csv", DataFrame(stack(dict2bifurcation(bfcn_lorenz63_rcvd)), ["h", "v"]))
bfcn_foodchain = JLD2.load("bifurcation_foodchain.jld2")["bfcn"]; CSV.write("bfcn.csv", DataFrame(stack(dict2bifurcation(bfcn_foodchain)), ["h", "v"]))
bfcn_foodchain_rcvd = JLD2.load("bifurcation_foodchain_recover.jld2")["bfcn"]; CSV.write("bfcn_.csv", DataFrame(stack(dict2bifurcation(bfcn_foodchain_rcvd)), ["h", "v"]))
bfcn_foodchain_fail = JLD2.load("bifurcation_foodchain_failed.jld2")["bfcn"]; CSV.write("bfcn__.csv", DataFrame(stack(dict2bifurcation(bfcn_foodchain_fail)), ["h", "v"]))
bfcn_thomas = JLD2.load("bifurcation_thomas.jld2")["bfcn"]; CSV.write("bfcn.csv", DataFrame(stack(dict2bifurcation(bfcn_thomas)), ["h", "v"]))
bfcn_thomas_rcvd = JLD2.load("bifurcation_thomas_recover.jld2")["bfcn"]; CSV.write("bfcn_.csv", DataFrame(stack(dict2bifurcation(bfcn_thomas_rcvd)), ["h", "v"]))
bfcn_thomas_fail = JLD2.load("bifurcation_thomas_failed.jld2")["bfcn"]; CSV.write("bfcn__.csv", DataFrame(stack(dict2bifurcation(bfcn_thomas_fail)), ["h", "v"]))
plot(
    scatter(dict2bifurcation(bfcn_lorenz63), color = :black, xticks = [0, 101, 1001], xformatter = _ -> ""; bfcnargs...),
    scatter(dict2bifurcation(bfcn_lorenz63_rcvd), color = :blue, xticks = [0, 1, 10]; bfcnargs...),
    layout = (:, 1), size= (400, 400)
); png("bifurcation_lorenz63.png")
plot(
    scatter(dict2bifurcation(bfcn_foodchain), color = :black, xticks = [.9, .95, .96, 1]; bfcnargs...),
    scatter(dict2bifurcation(bfcn_foodchain_rcvd), color = :blue, xticks = [-5, 0, 1, 5], xlims = [-5.3, 5.3]; bfcnargs...),
    scatter(dict2bifurcation(bfcn_foodchain_fail), color = :red, xticks = [-5, 0, 1, 5], xlims = [-5.3, 5.3]; bfcnargs...),
    layout = (:, 1), size = (400, 600)
); png("bifurcation_foodchain.png")
plot(
    scatter(dict2bifurcation(bfcn_thomas), color = :black, xticks = [.14, .141, .16], xformatter = _ -> ""; bfcnargs...),
    scatter(dict2bifurcation(bfcn_thomas_rcvd), color = :blue, xticks = [0, 1, 20]; bfcnargs...),
    scatter(dict2bifurcation(bfcn_thomas_fail), color = :red, xticks = [0, 1, 20]; bfcnargs...),
    layout = (:, 1), size = (400, 600)
); png("bifurcation_thomas.png")



bfcn1 = JLD2.load("bifurcation_foodchain_recover_950_955.jld2")["bfcn"]; CSV.write("bfcn1.csv", DataFrame(stack(dict2bifurcation(bfcn1)), ["h", "v"]))
bfcn2 = JLD2.load("bifurcation_foodchain_recover_950_970.jld2")["bfcn"]; CSV.write("bfcn2.csv", DataFrame(stack(dict2bifurcation(bfcn2)), ["h", "v"]))
bfcn3 = JLD2.load("bifurcation_foodchain_recover_930_970.jld2")["bfcn"]; CSV.write("bfcn3.csv", DataFrame(stack(dict2bifurcation(bfcn3)), ["h", "v"]))

