include.("../core/" .* readdir("core")[[1,2,3,4,6]])

using MATLAB

bfcnAh, bfcnAv = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnA.jld2"));

# mat"""
# plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
# """

bfcnB_9394h, bfcnB_9394v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnB_9394.jld2"));
bfcnB_9395h, bfcnB_9395v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnB_9395.jld2"));
bfcnB_9495h, bfcnB_9495v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnB_9495.jld2"));
bfcnB_9396h, bfcnB_9396v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnB_9396.jld2"));
bfcnB_9496h, bfcnB_9496v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnB_9496.jld2"));
bfcnB_9596h, bfcnB_9596v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnB_9596.jld2"));

bfcnC_9394h, bfcnC_9394v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnC_9394.jld2"));
bfcnC_9395h, bfcnC_9395v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnC_9395.jld2"));
bfcnC_9495h, bfcnC_9495v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnC_9495.jld2"));
bfcnC_9396h, bfcnC_9396v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnC_9396.jld2"));
bfcnC_9496h, bfcnC_9496v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnC_9496.jld2"));
bfcnC_9596h, bfcnC_9596v = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnC_9596.jld2"));

mat"""
nIDs = 26;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

fig = tiledlayout(4, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
set(gcf, 'Position', [100, 100, 900, 900])

axa1 = axes(fig);
axa1.Layout.Tile = 1;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('0.93   0.94    ', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axa1.XAxisLocation = 'top';
axa1.XLim = [0.88, 1.00]; axa1.YLim = [0.55, 0.8];
axa1.XTick = []; axa1.YTick = [0.55, 0.8];
axa2 = axes(fig);
axa2.Layout.Tile = 1;
plot($bfcnB_9394h, $bfcnB_9394v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
box on
grid on
axa2.XLim = [-5, 7]; axa2.YLim = [0.55, 0.8];
axa2.XTick = [-5, 0, 1, 7]; axa2.YTick = [];
axa2.Color = 'none';
text(0.02 , 1.05, charlbl{1}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axb1 = axes(fig);
axb1.Layout.Tile = 2;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('0.93     0.95', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axb1.XAxisLocation = 'top';
axb1.XLim = [0.88, 1.00]; axb1.YLim = [0.55, 0.8];
axb1.XTick = []; axb1.YTick = [0.55, 0.8];
axb2 = axes(fig);
axb2.Layout.Tile = 2;
plot($bfcnB_9395h, $bfcnB_9395v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
box on
grid on
axb2.XLim = [-2.5, 3.5]; axb2.YLim = [0.55, 0.8];
axb2.XTick = [-2.5, 0, 1, 3.5]; axb2.YTick = [];
axb2.Color = 'none';
text(0.02 , 1.05, charlbl{2}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axc1 = axes(fig);
axc1.Layout.Tile = 3;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('   0.94   0.95', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axc1.XAxisLocation = 'top';
axc1.XLim = [0.88, 1.00]; axc1.YLim = [0.55, 0.8];
axc1.XTick = []; axc1.YTick = [0.55, 0.8];
axc2 = axes(fig);
axc2.Layout.Tile = 3;
plot($bfcnB_9495h, $bfcnB_9495v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
box on
grid on
axc2.XLim = [-6, 6]; axc2.YLim = [0.55, 0.8];
axc2.XTick = [-6, 0, 1, 6]; axc2.YTick = [];
axc2.Color = 'none';
text(0.02 , 1.05, charlbl{3}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axd1 = axes(fig);
axd1.Layout.Tile = 4;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('   0.93        0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axd1.XAxisLocation = 'top';
axd1.XLim = [0.88, 1.00]; axd1.YLim = [0.55, 0.8];
axd1.XTick = []; axd1.YTick = [0.55, 0.8];
axd2 = axes(fig);
axd2.Layout.Tile = 4;
plot($bfcnB_9396h, $bfcnB_9396v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
box on
grid on
axd2.XLim = [-1.66, 2.33]; axd2.YLim = [0.55, 0.8];
axd2.XTick = [-1.66, 0, 1, 2.33]; axd2.YTick = [];
axd2.XTickLabel = [-1.6, 0, 1, 2.3];
axd2.Color = 'none';
text(0.02 , 1.05, charlbl{4}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axe1 = axes(fig);
axe1.Layout.Tile = 5;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('         0.94   0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axe1.XAxisLocation = 'top';
axe1.XLim = [0.88, 1.00]; axe1.YLim = [0.55, 0.8];
axe1.XTick = []; axe1.YTick = [0.55, 0.8];
axe2 = axes(fig);
axe2.Layout.Tile = 5;
plot($bfcnB_9496h, $bfcnB_9496v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
box on
grid on
axe2.XLim = [-3, 3]; axe2.YLim = [0.55, 0.8];
axe2.XTick = [-3, 0, 1, 3]; axe2.YTick = [];
axe2.Color = 'none';
text(0.02 , 1.05, charlbl{5}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axf1 = axes(fig);
axf1.Layout.Tile = 6;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('               0.95   0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axf1.XAxisLocation = 'top';
axf1.XLim = [0.88, 1.00]; axf1.YLim = [0.55, 0.8];
axf1.XTick = []; axf1.YTick = [0.55, 0.8];
axf2 = axes(fig);
axf2.Layout.Tile = 6;
plot($bfcnB_9596h, $bfcnB_9596v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red')
box on
grid on
axf2.XLim = [-7, 5]; axf2.YLim = [0.55, 0.8];
axf2.XTick = [-7, 0, 1, 5]; axf2.YTick = [];
axf2.Color = 'none';
text(0.02 , 1.05, charlbl{6}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axg1 = axes(fig);
axg1.Layout.Tile = 7;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('0.93   0.94    ', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axg1.XAxisLocation = 'top';
axg1.XLim = [0.88, 1.00]; axg1.YLim = [0.55, 0.8];
axg1.XTick = []; axg1.YTick = [0.55, 0.8];
axg2 = axes(fig);
axg2.Layout.Tile = 7;
plot($bfcnC_9394h, $bfcnC_9394v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
box on
grid on
axg2.XLim = [-5, 7]; axg2.YLim = [0.55, 0.8];
axg2.XTick = [-5, 0, 1, 7]; axg2.YTick = [];
axg2.Color = 'none';
text(0.02 , 1.05, charlbl{7}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axh1 = axes(fig);
axh1.Layout.Tile = 8;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('0.93     0.95', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axh1.XAxisLocation = 'top';
axh1.XLim = [0.88, 1.00]; axh1.YLim = [0.55, 0.8];
axh1.XTick = []; axh1.YTick = [0.55, 0.8];
axh2 = axes(fig);
axh2.Layout.Tile = 8;
plot($bfcnC_9395h, $bfcnC_9395v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
box on
grid on
axh2.XLim = [-2.5, 3.5]; axh2.YLim = [0.55, 0.8];
axh2.XTick = [-2.5, 0, 1, 3.5]; axh2.YTick = [];
axh2.Color = 'none';
text(0.02 , 1.05, charlbl{8}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axi1 = axes(fig);
axi1.Layout.Tile = 9;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('   0.94   0.95', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axi1.XAxisLocation = 'top';
axi1.XLim = [0.88, 1.00]; axi1.YLim = [0.55, 0.8];
axi1.XTick = []; axi1.YTick = [0.55, 0.8];
axi2 = axes(fig);
axi2.Layout.Tile = 9;
plot($bfcnC_9495h, $bfcnC_9495v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
box on
grid on
axi2.XLim = [-6, 6]; axi2.YLim = [0.55, 0.8];
axi2.XTick = [-6, 0, 1, 6]; axi2.YTick = [];
axi2.Color = 'none';
text(0.02 , 1.05, charlbl{9}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axj1 = axes(fig);
axj1.Layout.Tile = 10;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('   0.93        0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axj1.XAxisLocation = 'top';
axj1.XLim = [0.88, 1.00]; axj1.YLim = [0.55, 0.8];
axj1.XTick = []; axj1.YTick = [0.55, 0.8];
axj2 = axes(fig);
axj2.Layout.Tile = 10;
plot($bfcnC_9396h, $bfcnC_9396v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
box on
grid on
axj2.XLim = [-1.66, 2.33]; axj2.YLim = [0.55, 0.8];
axj2.XTick = [-1.66, 0, 1, 2.33]; axj2.YTick = [];
axj2.XTickLabel = [-1.6, 0, 1, 2.3];
axj2.Color = 'none';
text(0.02 , 1.05, charlbl{10}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axk1 = axes(fig);
axk1.Layout.Tile = 11;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('         0.94   0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axk1.XAxisLocation = 'top';
axk1.XLim = [0.88, 1.00]; axk1.YLim = [0.55, 0.8];
axk1.XTick = []; axk1.YTick = [0.55, 0.8];
axk2 = axes(fig);
axk2.Layout.Tile = 11;
plot($bfcnC_9496h, $bfcnC_9496v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
box on
grid on
axk2.XLim = [-3, 3]; axk2.YLim = [0.55, 0.8];
axk2.XTick = [-3, 0, 1, 3]; axk2.YTick = [];
axk2.Color = 'none';
text(0.02 , 1.05, charlbl{11}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axl1 = axes(fig);
axl1.Layout.Tile = 12;
plot($bfcnAh, $bfcnAv, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
title('               0.95   0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axl1.XAxisLocation = 'top';
axl1.XLim = [0.88, 1.00]; axl1.YLim = [0.55, 0.8];
axl1.XTick = []; axl1.YTick = [0.55, 0.8];
axl2 = axes(fig);
axl2.Layout.Tile = 12;
plot($bfcnC_9596h, $bfcnC_9596v, 'o', 'MarkerSize', 0.5, 'LineStyle', 'none', 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue')
box on
grid on
axl2.XLim = [-7, 5]; axl2.YLim = [0.55, 0.8];
axl2.XTick = [-7, 0, 1, 5]; axl2.YTick = [];
axl2.Color = 'none';
text(0.02 , 1.05, charlbl{12}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
"""

bfcnAh, bfcnAv = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnA.jld2"));

function tickfmt(x)
    if isapprox(x, round(x); atol=1e-8)
        return string(Int(round(x)))
    else
        s = string(round(x, digits=1))
        return s
    end
end

β0, β1 = 0, 1
pm, pM = 0.88, 1.0

pp_ = [[93, 94], [93, 95], [94, 95], [93, 96], [94, 96], [95, 96]]

plt_ = []
for (p0, p1) = pp_
bfcnh, bfcnv = dict2bifurcation(callbfcn("G:/BF/foodchain/-4/bfcnB_$(p0)$(p1).jld2"))
p0 /= 100; p1 /= 100;
βm = (pm - p0) / (p1 - p0)
βM = (pM - p0) / (p1 - p0)
plt = scatter(bfcnAh, bfcnAv, color = :gray, msw = 0, ms = 0.3, ylims = [0.55, 0.8], yticks = ([0.55, 0.8], ["0.55", "0.8"]), xticks = ([pm, p0, p1, pM], tickfmt.([βm, β0, β1, βM])), xlims = [pm, pM]);
annotate!(plt, (-0.1, 0.5), text(L"P", 12, rotation = 90))
annotate!(plt, (0.25, -0.1), text(L"\beta", 12))
scatter!(twiny(plt), bfcnh, bfcnv, color = :red, msw = 0, ms = 0.3, ylims = [0.55, 0.8], xticks = ([β0, β1], ["$p0   ", "   $p1"]), xlims = [βm, βM], grid = true);
scatter!(twinx(plt), yticks = []); 
push!(plt_, plt)
end

plot(plt_..., layout = (:, 3), size = [900, 450]); @time png("temp4")
