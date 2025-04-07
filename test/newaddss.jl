include("../core/header.jl")


##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
# cnfg = (; f_ = [cospi], λ = 1e-1)
dt = 1e-5; tspan = [0, 10]; θ = 1e-16;
# dt = 1e-5 is good to bifurcation diagram but not for lyapunov spectrum, 1e-6 is needed
# θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
if isfile("hyperparameter/0 data.csv")
    data = CSV.read("hyperparameter/0 data.csv", DataFrame)
    # add_subsystem!(data, vrbl, (λ = 1e-1, N = 1, M = 1); θ)
else
    @time data = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
    # add_subsystem!(data, vrbl, (λ = 1e-1, N = 1, M = 1); θ)
    CSV.write("hyperparameter/0 data.csv", data, bom = true)
end
# plot(data.u, color = data.subsystem)


@show jumpt = detect_jump(data, vrbl);

sets = set_divider(jumpt)
cnfg = (; N = 2, M = 2)
datasets = [data[set, :] for set in sets]
# f_ = [SINDy(dataset, vrbl...; cnfg...) for dataset in datasets]
# MSE_ = [f_.MSE for f_ in f_]
label = zeros(Int64, nrow(data));
subsystem = zeros(Int64, nrow(data));
idcs = eachindex(sets)
volume = length.(sets)

# ----------------------------------------------

idx_hero = argmax(volume)
aliens = setdiff(idcs, idx_hero .+ [-1, 0, 1])

report_hero = DataFrame(alien = Int64[], f = [], mse = Float64[])
for alien = aliens
    data_alien = [datasets[[alien, idx_hero]]...;]
    f = SINDy(data_alien, vrbl...; cnfg...)
    push!(report_hero, [alien, f, f.MSE])
end
scatter(report_hero.alien, report_hero.mse, yscale = :log10)
worst = SINDy([datasets[idcs]...;], vrbl...; cnfg...)
worst.MSE

idx_hive = report_hero.alien[argmax(report_hero.mse)]
report_hive = DataFrame(alien = Int64[], f = [], mse = Float64[])
for alien = aliens
    data_hive = [datasets[[alien, idx_hive]]...;]
    f = SINDy(data_hive, vrbl...; cnfg...)
    push!(report_hive, [alien, f, f.MSE])
end
# scatter(report_hive.alien, report_hive.mse, yscale = :log10)


tempidx = report_hero.mse .< report_hive.mse
scatter(report_hero.alien, report_hero.mse, yscale = :log10, xticks = idcs)
scatter!(report_hero.alien[tempidx], report_hero.mse[tempidx], yscale = :log10, color = :red, shape = :x)

# ----------------------------------------------
idx_labeled = [idx_hero; report_hero.alien[tempidx]]
volume[idx_labeled] .= 0
idcs = setdiff(idcs, idx_labeled)

idx_hero = argmax(volume)
aliens = setdiff(idcs, idx_hero .+ [-1, 0, 1])

report_hero = DataFrame(alien = Int64[], f = [], mse = Float64[])
for alien = aliens
    data_alien = [datasets[[alien, idx_hero]]...;]
    f = SINDy(data_alien, vrbl...; cnfg...)
    push!(report_hero, [alien, f, f.MSE])
end
scatter(report_hero.alien, report_hero.mse, yscale = :log10)
worst = SINDy([datasets[idcs]...;], vrbl...; cnfg...)
worst.MSE

idx_hive = report_hero.alien[argmax(report_hero.mse)]
report_hive = DataFrame(alien = Int64[], f = [], mse = Float64[])
for alien = aliens
    data_hive = [datasets[[alien, idx_hive]]...;]
    f = SINDy(data_hive, vrbl...; cnfg...)
    push!(report_hive, [alien, f, f.MSE])
end
# scatter(report_hive.alien, report_hive.mse, yscale = :log10)


tempidx = report_hero.mse .< report_hive.mse
scatter(report_hero.alien, report_hero.mse, yscale = :log10, xticks = idcs)
scatter!(report_hero.alien[tempidx], report_hero.mse[tempidx], yscale = :log10, color = :red, shape = :x)


# ----------------------------------------------
idx_labeled = [idx_hero; report_hero.alien[tempidx]]
volume[idx_labeled] .= 0
idcs = setdiff(idcs, idx_labeled)

idx_hero = argmax(volume)
aliens = setdiff(idcs, idx_hero .+ [-1, 0, 1])

report_hero = DataFrame(alien = Int64[], f = [], mse = Float64[])
for alien = aliens
    data_alien = [datasets[[alien, idx_hero]]...;]
    f = SINDy(data_alien, vrbl...; cnfg...)
    push!(report_hero, [alien, f, f.MSE])
end
scatter(report_hero.alien, report_hero.mse, yscale = :log10)
worst = SINDy([datasets[idcs]...;], vrbl...; cnfg...)
worst.MSE

report_hero.mse

# -----------------------------------------------

