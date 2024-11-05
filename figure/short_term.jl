include("../core/header.jl")


##########################################################################
#                                                                        #
#                            Soft impact model                           #
#                                                                        #
##########################################################################
vrbl = [:dt, :du, :dv], [:t, :u, :v]
cnfg = (; f_ = [cospi, sign], λ = 1e-1)
dt = 1e-5; θ = 1e-16;
# dt = 1e-5 is good to bifurcation diagram but not for lyapunov spectrum, 1e-6 is needed
# θ1 = 1e-8; θ2 = 1e-12; θ3 = 1e-5; min_rank = 21;
@time trng = factory_soft(DataFrame, 0.1, tspan = [0, 10], ic = [0, .0446272, -0.119564]; dt)
add_subsystem!(trng, vrbl, cnfg; θ)
f_ = [SINDy(df, vrbl...; cnfg...) for df in groupby(trng, :subsystem)] # print.(f_)
Dtree = dryad(trng, last(vrbl)); # print_tree(Dtree)
prd1 = DataFrame(solve(f_, collect(trng[1, last(vrbl)]), dt, 0:dt:10, Dtree), last(vrbl))

@time test = factory_soft(DataFrame, 0.1, tspan = [0, 10], ic = [0, .000129715, .301469]; dt)
prd2 = DataFrame(solve(f_, collect(test[1, last(vrbl)]), dt, 0:dt:10, Dtree), last(vrbl))

q1 = plot(xticks = [0, 10], yticks = [.05, -.05], xlims = [0, 10], legend = :none)
plot!(trng.t[1:1000:end], trng.u[1:1000:end], alpha = .5, lw = 2, color = :black)
plot!(prd1.t[1:1000:end], prd1.u[1:1000:end], alpha = .5, lw = 2, color = 1)

q2 = plot(xticks = [0, 10], yticks = [], xlims = [0, 10], legend = :none)
plot!(trng.t[1:1000:end], trng.v[1:1000:end], alpha = .5, lw = 2, color = :black)
plot!(prd1.t[1:1000:end], prd1.v[1:1000:end], alpha = .5, lw = 2, color = 1)

q3 = plot(xticks = [0, 10], yticks = [.05, -.05], xlims = [0, 10], legend = :none)
plot!(test.t[1:1000:end], test.u[1:1000:end], alpha = .5, lw = 2, color = :black)
plot!(prd2.t[1:1000:end], prd2.u[1:1000:end], alpha = .5, lw = 2, color = 1)

q4 = plot(xticks = [0, 10], yticks = [], xlims = [0, 10], legend = :none)
plot!(test.t[1:1000:end], test.v[1:1000:end], alpha = .5, lw = 2, color = :black)
plot!(prd2.t[1:1000:end], prd2.v[1:1000:end], alpha = .5, lw = 2, color = 1)

plot(q1, q2, q3, q4, layout = (4, 1), size = (800, 800), right_margin = 3mm, dpi = 300)
png("q1234")

@info "------------------------------------------------------------------"

# data_buck = CSV.read("cached_buck.csv", DataFrame)[1:50000,:]
# blt_b0 = plot(xticks = [0, 0.005])
# plt_b1 = plot(blt_b0, data_buck.t, data_buck.V, yticks = [11.75, 12.75])
# plt_b1 = plot!(data_buck.t, data_buck.Vr, color = 1)
# plt_b2 = plot(blt_b0, data_buck.t, data_buck.I, yticks = [.4, .7])
# plt_b3 = plot(data_buck.V, data_buck.I, xticks = [11.75, 12.75], yticks = [.4, .7])
# lyt_b0 = @layout [[a; b] c]
# plot(plt_b1, plt_b2, plt_b3, layout = lyt_b0, size = [900, 500])
# png("eda_buck")

# data_soft = CSV.read("cached_soft.csv", DataFrame)[1:10:500000,:]
# blt_s0 = plot(xticks = [0, 5])
# plt_s1 = plot(blt_s0, data_soft.t, data_soft.u, yticks = [-.05, .05])
# plt_s1 = hline!([-.05, .05], color = 1)
# plt_s2 = plot(blt_s0, data_soft.t, data_soft.v, yticks = [-.4, .4])
# plt_s3 = plot(data_soft.u, data_soft.v)
# plt_s3 = vline!([-.05, .05], color = 1, xticks = [-.05, .05], yticks = [-.4, .4])
# lyt_s0 = @layout [[a; b] c]
# plot(plt_s1, plt_s2, plt_s3, layout = lyt_s0, size = [900, 500])
# png("eda_soft")

# data_hrnm = CSV.read("cached_hrnm.csv", DataFrame)[1:10:100000,:]
# blt_h0 = plot(xticks = [0, 100])
# plt_h1 = plot(blt_h0, data_hrnm.t, data_hrnm.x, yticks = [-1, 2])
# plt_h2 = plot(blt_h0, data_hrnm.t, data_hrnm.y, yticks = [-12, 0])
# plt_h3 = plot(blt_h0, data_hrnm.t, data_hrnm.z, yticks = [-4, 2])
# plt_h3 = hline!([-1, 1], color = 1)
# plt_h4 = plot(data_hrnm.x, data_hrnm.y, data_hrnm.z, xticks = [-1, 2], yticks = [-12, 0], zticks = [-4, 2])
# lyt_h0 = @layout [[a; b; d] c]
# plot(plt_h1, plt_h2, plt_h3, plt_h4, layout = lyt_h0, size = [900, 500])
# png("eda_hrnm")


using L1TrendFiltering
y = trng.v[1:100:end]
@time x = l1tf(y, 1).x;
plot(y, label = "ground truth"); plot!(x, label = "l1 trend filtering")
plot(abs.(diff(x)), label = "∇x, l1 trend filtering")
plot(diff(y), label = "∇x, l1 trend filtering")

threshold = 1e-4; 
D = diff(speye(length(x)), 2);
BP_L1 = find(round(D * y_pred_L1, 1))+1;
BP_L1 = countBP(BP_L1);

function bm(data, vrbl, cnfg;)
    normeddf = norm.(eachrow(diff(diff(Matrix(data[:, first(vrbl)]), dims = 1), dims = 1))) # scatter(normeddf[1:100:end], yscale = :log10)
    _jumpt = [0]; jumpt = [0];
    while true
        idx = argmax(normeddf)
        if all(abs.(_jumpt .- idx) .> 1)
            jumpt = deepcopy(_jumpt)
            push!(_jumpt, idx, idx+1, idx-1)
            normeddf[[idx, idx+1, idx-1]] .= -Inf
        else
            break
        end
    end
    jumpt = [1; sort(jumpt[2:end]); nrow(data)]
    sets = set_divider(jumpt)
    return jumpt
end
@time result = bm(trng[1:100:end,:], vrbl, cnfg);


plot(y, label = "ground truth"); scatter!(result, trng.v[result], label = "proposed")