include("../core/visual.jl")

data_buck = CSV.read("cached_buck.csv", DataFrame)[1:50000,:]
blt_b0 = plot(xticks = [0, 0.005])
plt_b1 = plot(blt_b0, data_buck.t, data_buck.V, yticks = [11.75, 12.75])
plt_b1 = plot!(data_buck.t, data_buck.Vr, color = :blue)
plt_b2 = plot(blt_b0, data_buck.t, data_buck.I, yticks = [.4, .7])
plt_b3 = plot(data_buck.V, data_buck.I, xticks = [11.75, 12.75], yticks = [.4, .7])
lyt_b0 = @layout [[a; b] c]
plot(plt_b1, plt_b2, plt_b3, layout = lyt_b0, size = [900, 500])
png("eda_buck")

data_soft = CSV.read("cached_soft.csv", DataFrame)[1:10:500000,:]
blt_s0 = plot(xticks = [0, 5])
plt_s1 = plot(blt_s0, data_soft.t, data_soft.u, yticks = [-.05, .05])
plt_s1 = hline!([-.05, .05], color = :blue)
plt_s2 = plot(blt_s0, data_soft.t, data_soft.v, yticks = [-.4, .4])
plt_s3 = plot(data_soft.u, data_soft.v)
plt_s3 = vline!([-.05, .05], color = :blue, xticks = [-.05, .05], yticks = [-.4, .4])
lyt_s0 = @layout [[a; b] c]
plot(plt_s1, plt_s2, plt_s3, layout = lyt_s0, size = [900, 500])
png("eda_soft")

data_hrnm = CSV.read("cached_hrnm.csv", DataFrame)[1:10:100000,:]
blt_h0 = plot(xticks = [0, 100])
plt_h1 = plot(blt_h0, data_hrnm.t, data_hrnm.x, yticks = [-1, 2])
plt_h2 = plot(blt_h0, data_hrnm.t, data_hrnm.y, yticks = [-12, 0])
plt_h3 = plot(blt_h0, data_hrnm.t, data_hrnm.z, yticks = [-4, 2])
plt_h3 = hline!([-1, 1], color = :blue)
plt_h4 = plot(data_hrnm.x, data_hrnm.y, data_hrnm.z, xticks = [-1, 2], yticks = [-12, 0], zticks = [-4, 2])
lyt_h0 = @layout [[a; b; d] c]
plot(plt_h1, plt_h2, plt_h3, plt_h4, layout = lyt_h0, size = [900, 500])
png("eda_hrnm")
