include("../core/header.jl")


##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 1e-1)
dt = 1e-5; tspan = [0, 10]; θ = 1e-16;
# dt = 1e-5 is good to bifurcation diagram but not for lyapunov spectrum, 1e-6 is needed
# θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
@time trng = factory_soft(DataFrame, 0.1, ic = [0, .0446272, -0.119564]; tspan, dt)
add_subsystem!(trng, vrbl, cnfg; θ)
f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(trng, :subsystem)] # print.(f_)
Dtree = dryad(trng, last(vrbl)); # print_tree(Dtree)
# prd1 = DataFrame(solve(f_, collect(trng[1, last(vrbl)]), dt, 0:dt:10, Dtree), last(vrbl))

@time test = factory_soft(DataFrame, 0.1, ic = [0, .000129715, .301469]; tspan = 2tspan, dt)
prd2 = DataFrame(solve(f_, collect(test[1, last(vrbl)]), dt, first(2tspan):dt:last(2tspan), Dtree), last(vrbl))[1:end-1, :]

# q1 = plot(xticks = [0, 10], yticks = [.05, -.05], xlims = [0, 10], legend = :none)
# plot!(trng.t[1:1000:end], trng.u[1:1000:end], alpha = .5, lw = 2, color = :black)
# plot!(prd1.t[1:1000:end], prd1.u[1:1000:end], alpha = .5, lw = 2, color = 1)

# q2 = plot(xticks = [0, 10], yticks = [], xlims = [0, 10], legend = :none)
# plot!(trng.t[1:1000:end], trng.v[1:1000:end], alpha = .5, lw = 2, color = :black)
# plot!(prd1.t[1:1000:end], prd1.v[1:1000:end], alpha = .5, lw = 2, color = 1)

# q3 = plot(xticks = [0, 10], yticks = [.05, -.05], xlims = [0, 10], legend = :none)
# plot!(test.t[1:1000:end], test.u[1:1000:end], alpha = .5, lw = 2, color = :black)
# plot!(prd2.t[1:1000:end], prd2.u[1:1000:end], alpha = .5, lw = 2, color = 1)

# q4 = plot(xticks = [0, 10], yticks = [], xlims = [0, 10], legend = :none)
# plot!(test.t[1:1000:end], test.v[1:1000:end], alpha = .5, lw = 2, color = :black)
# plot!(prd2.t[1:1000:end], prd2.v[1:1000:end], alpha = .5, lw = 2, color = 1)

# plot(q1, q2, q3, q4, layout = (4, 1), size = (800, 800), right_margin = 3mm, dpi = 300)
# png("q1234")
@time temp1 = factory_soft(DataFrame, 0.1, ic = [0, .000129715, .301469]; tspan = 10tspan, dt)
temp2 = DataFrame(solve(f_, collect(temp1[1, last(vrbl)]), dt, first(10tspan):dt:last(10tspan), Dtree), last(vrbl))[1:end-1, :]

uvargs = (; xlabel = L"u", ylabel = L"\dot{u}", legend = :none, lw = .1)
plot(
    plot(temp1.u[1:100:end], temp1.v[1:100:end], title = "reference", color = :blue; uvargs...),
    plot(temp2.u[1:100:end], temp2.v[1:100:end], title = "recovered", color = :red ; uvargs...),
    size = [600, 300], dpi = 300
    ); png("uv")
@info "------------------------------------------------------------------"

rename!(test, "t_".* names(test))
rename!(prd2, "p_".* names(prd2))
CSV.write("soft_shortterm.csv", [test prd2][1:1000:end, :])

##########################################################################
#                                                                        #
#                             Gear system                                #
#                                                                        #
##########################################################################
vrbl = [:dx, :dv, :dΩ, :dθ], [:x, :v, :Ω, :θ]
cnfg = (; N = 1, f_ = [cos], C = 2,  λ = 1e-2)
dt = 1e-2; tspan = [0, 100]; θ = 1e-20; dos = 1

@time trng = factory_gear(DataFrame, -0.11; tspan, dt)
add_subsystem!(trng, vrbl, cnfg; θ, dos)
f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(trng, :subsystem)] # print.(f_)
Dtree = dryad(trng, last(vrbl)); # print_tree(Dtree)
# prd1 = DataFrame(solve(f_, collect(trng[1, last(vrbl)]), dt, 0:dt:10, Dtree), last(vrbl))

@time test = factory_gear(DataFrame, -0.11, ic = [0.2, 0.2, 0.1, 0.0]; tspan = 3tspan, dt)
prd2 = DataFrame(solve(f_, collect(test[1, last(vrbl)]), dt, first(3tspan):dt:last(3tspan), Dtree), last(vrbl))[1:end-1, :]

rename!(test, "t_".* names(test))
rename!(prd2, "p_".* names(prd2))
CSV.write("gear_shortterm.csv", [test prd2], bom = true)


##########################################################################
#                                                                        #
#                           Hindmarsh-Rose model                         #
#                                                                        #
##########################################################################
vrbl = [:dt, :dx, :dy, :dz], [:t, :x, :y, :z]
cnfg = (; N = 3, f_ = [cos])
dt = 1e-3; θ = 1e-20;

@time trng = factory_hrnm(DataFrame, 0.1; tspan, dt)
add_subsystem!(trng, vrbl, cnfg; θ)
f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(trng, :subsystem)] # print.(f_)
Dtree = dryad(trng, last(vrbl)); # print_tree(Dtree)
# prd1 = DataFrame(solve(f_, collect(trng[1, last(vrbl)]), dt, 0:dt:10, Dtree), last(vrbl))

@time test = factory_hrnm(DataFrame, 0.1, ic = [0, 0.1, 0.2, 0.1]; tspan = 2tspan, dt)
prd2 = DataFrame(solve(f_, collect(test[1, last(vrbl)]), dt, first(2tspan):dt:last(2tspan), Dtree), last(vrbl))[1:end-1, :]

# plot(test.x, color = :black)
# plot!(prd2.x, color = :blue)

rename!(test, "t_".* names(test))
rename!(prd2, "p_".* names(prd2))
CSV.write("hrnm_shortterm.csv", [test prd2], bom = true)


plot(short.t, short.z, color = :black)
short = data[173000:10:195000, :]
pargs = (; ticks = false, legend = :none, dpi = 300, size = [800, 200], lw = 2, ylims = [-2.5, 1.7])
plot(short.t, short.z, color = :black; pargs...); png("1")
plot(short.t, short.z, color = short.subsystem; pargs...); png("2")
