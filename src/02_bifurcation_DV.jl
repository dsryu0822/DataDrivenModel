include.("../core/" .* readdir("core")[[1,2,3,4,6]])

function box_count(V, ticks)
    bin_ = []
    for j in axes(V, 2)
        tick = ticks[j]
        v_ = V[:, j]
        push!(bin_, sum([(v_ .≤ t) for t in tick]))
    end
    freq = zeros(Int64, (length.(ticks))...)
    for ijk in zip(bin_...)
        @inbounds freq[ijk...] += 1
    end
    return freq
end
# box_count(randn(1000, 3), [-2:1e-2:2, -2:1e-2:2, -2:1e-2:2])
# V = traj1
# ticks = tick_

function dv(traj1, traj2)
    if isempty(traj1) || isempty(traj2)
        return NaN
    end
    if !all(isfinite.(Matrix(traj2)))
        return 2.0
    end
    try
        d = size(traj1, 2)
        min_ = minimum.(eachcol(vcat(traj1, traj2)))
        max_ = maximum.(eachcol(vcat(traj1, traj2)))
        tick_ = [range(min_[i] - eps(), max_[i] + eps(), length = 11) for i in 1:d]
        freq1 = box_count(traj1, tick_)
        freq2 = box_count(traj2, tick_)
        freq1 = freq1 / sum(freq1)
        freq2 = freq2 / sum(freq2)
        return sum(abs, freq1 - freq2)        
    catch e
        return NaN
    end
end


resultDV = DataFrame(zeros(2001, 13), ["k", "B9394", "B9395", "B9495", "B9396", "B9496", "B9596", "C9394", "C9395", "C9495", "C9396", "C9496", "C9596"])
@showprogress @threads for k = 1:2001
    idx = lpad(k, 4, '0')
    trajA =      try CSV.read("G:/BF/foodchain/trajA/$(idx).csv", DataFrame)[:, [:R, :C, :P]]      catch e DataFrame() end
    trajB_9394 = try CSV.read("G:/BF/foodchain/trajB_9394/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajB_9395 = try CSV.read("G:/BF/foodchain/trajB_9395/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajB_9495 = try CSV.read("G:/BF/foodchain/trajB_9495/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajB_9396 = try CSV.read("G:/BF/foodchain/trajB_9396/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajB_9496 = try CSV.read("G:/BF/foodchain/trajB_9496/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajB_9596 = try CSV.read("G:/BF/foodchain/trajB_9596/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajC_9394 = try CSV.read("G:/BF/foodchain/trajC_9394/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajC_9395 = try CSV.read("G:/BF/foodchain/trajC_9395/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajC_9495 = try CSV.read("G:/BF/foodchain/trajC_9495/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajC_9396 = try CSV.read("G:/BF/foodchain/trajC_9396/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajC_9496 = try CSV.read("G:/BF/foodchain/trajC_9496/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    trajC_9596 = try CSV.read("G:/BF/foodchain/trajC_9596/$(idx).csv", DataFrame)[:, [:R, :C, :P]] catch e DataFrame() end
    resultDV[k, :] = [k, dv(trajA, trajB_9394), dv(trajA, trajB_9395), dv(trajA, trajB_9495), dv(trajA, trajB_9396), dv(trajA, trajB_9496), dv(trajA, trajB_9596), dv(trajA, trajC_9394), dv(trajA, trajC_9395), dv(trajA, trajC_9495), dv(trajA, trajC_9396), dv(trajA, trajC_9496), dv(trajA, trajC_9596)]
end
CSV.write("G:/BF/foodchain/DV.csv", resultDV)

plt_ = []
for col in eachcol(resultDV)[2:end]
    push!(plt_, plot(col))
end
plot(plt_..., ylims = [0, 0.05], layout = (4, 3), size = (900, 900))


using MATLAB

mat"""
tiledlayout(4, 3, 'Padding', 'compact');

nexttile
plot(linspace(-5, 7, 2001), $(resultDV.B9394), 'r')
xlim([-5, 7]); xticks([-5, 0, 1, 7]); ylim([0, 2]);

nexttile
plot(linspace(-2.5, 3.5, 2001), $(resultDV.B9395), 'r')
xlim([-2.5, 3.5]); xticks([-2.5, 0, 1, 3.5]); ylim([0, 2]);

nexttile
plot(linspace(-6, 6, 2001), $(resultDV.B9495), 'r')
xlim([-6, 6]); xticks([-6, 0, 1, 6]); ylim([0, 2]);

nexttile
plot(linspace(-1.7, 2.3, 2001), $(resultDV.B9396), 'r')
xlim([-1.7, 2.3]); xticks([-1.7, 0, 1, 2.3]); ylim([0, 2]);

nexttile
plot(linspace(-3, 3, 2001), $(resultDV.B9496), 'r')
xlim([-3, 3]); xticks([-3, 0, 1, 3]); ylim([0, 2]);

nexttile
plot(linspace(-7, 5, 2001), $(resultDV.B9596), 'r')
xlim([-7, 5]); xticks([-7, 0, 1, 5]); ylim([0, 2]);

nexttile
plot(linspace(-5, 7, 2001), $(resultDV.C9394), 'b')
xlim([-5, 7]); xticks([-5, 0, 1, 7]); ylim([0, 2]);

nexttile
plot(linspace(-2.5, 3.5, 2001), $(resultDV.C9395), 'b')
xlim([-2.5, 3.5]); xticks([-2.5, 0, 1, 3.5]); ylim([0, 2]);

nexttile
plot(linspace(-6, 6, 2001), $(resultDV.C9495), 'b')
xlim([-6, 6]); xticks([-6, 0, 1, 6]); ylim([0, 2]);

nexttile
plot(linspace(-1.7, 2.3, 2001), $(resultDV.C9396), 'b')
xlim([-1.7, 2.3]); xticks([-1.7, 0, 1, 2.3]); ylim([0, 2]);

nexttile
plot(linspace(-3, 3, 2001), $(resultDV.C9496), 'b')
xlim([-3, 3]); xticks([-3, 0, 1, 3]); ylim([0, 2]);

nexttile
plot(linspace(-7, 5, 2001), $(resultDV.C9596), 'b')
xlim([-7, 5]); xticks([-7, 0, 1, 5]); ylim([0, 2]);
"""
