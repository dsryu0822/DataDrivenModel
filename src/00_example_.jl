include.("../core/" .* readdir("core")[[1,2,3,4,6]])

using MATLAB

bfcnAAh, bfcnAAv = dict2bifurcation(callbfcn("G:/BF/aizawa/bfcnA.jld2"));
bfcnBAh, bfcnBAv = dict2bifurcation(callbfcn("G:/BF/aizawa/bfcnB.jld2"));

bfcnABh, bfcnABv = dict2bifurcation(callbfcn("G:/BF/bouali/bfcnA.jld2"));
bfcnBBh, bfcnBBv = dict2bifurcation(callbfcn("G:/BF/bouali/bfcnB.jld2"));

bfcnADh, bfcnADv = dict2bifurcation(callbfcn("G:/BF/dadras/bfcnA.jld2"));
bfcnBDh, bfcnBDv = dict2bifurcation(callbfcn("G:/BF/dadras/bfcnB.jld2"));

bfcnAFh, bfcnAFv = dict2bifurcation(callbfcn("G:/BF/fourwing/bfcnA.jld2"));
bfcnBFh, bfcnBFv = dict2bifurcation(callbfcn("G:/BF/fourwing/bfcnB.jld2"));

bfcnARh, bfcnARv = dict2bifurcation(callbfcn("G:/BF/rikitake/bfcnA.jld2"));
bfcnBRh, bfcnBRv = dict2bifurcation(callbfcn("G:/BF/rikitake/bfcnB.jld2"));

mat"""
nIDs = 26;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

fig = tiledlayout(5, 3, 'Padding', 'compact');
set(gcf, 'Position', [100, 100, 900, 900])


nexttile
plot($bfcnAAh, $bfcnAAv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
xlim([0.6, 0.9]); ylim([0.5, 2]); xticks([0.6, 0.72, 0.76, 0.9]); yticks([0.5, 2]);
grid on
text(0, 1.06, charlbl{1}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot($bfcnBAh, $bfcnBAv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
xlim([-3, 4.5]); ylim([0.5, 2]); xticks([-3, 0, 1, 4.5]); yticks([0.5, 2]);
grid on

axa1 = axes(fig);
axa1.Layout.Tile = 3;
plot($bfcnAAh, $bfcnAAv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
xticks([]); yticks([]); xlim([0.6, 0.9]); ylim([0.5, 2]);
hold on
axa2 = axes(fig);
axa2.Layout.Tile = 3;
plot($bfcnBAh, $bfcnBAv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
xticks([-3, 0, 1, 4.5]); yticks([0.5, 2]); xlim([-3, 4.5]); ylim([0.5, 2]);
axa2.Color = 'none';
grid on



pm = 0.5; pM = 1.1; p0 = 0.7; p1 = 0.8;
bm = (pm - p0) / (p1 - p0);
bM = (pM - p0) / (p1 - p0);

nexttile
plot($bfcnABh, $bfcnABv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
xlim([pm, pM]); ylim([-2, 6]); xticks([pm, p0, p1, pM]); yticks([-2, 6]);
grid on
text(0, 1.06, charlbl{2}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot($bfcnBBh, $bfcnBBv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
xlim([bm, bM]); ylim([-2, 6]); xticks([bm, 0, 1, bM]); yticks([-2, 6]);
grid on

axb1 = axes(fig);
axb1.Layout.Tile = 6;
plot($bfcnABh, $bfcnABv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
axb1.XLim = [pm, pM]; axb1.YLim = [-2, 6];
axb1.XTick = []; axb1.YTick = [-2, 6];
hold on
axb2 = axes(fig);
axb2.Layout.Tile = 6;
plot($bfcnBBh, $bfcnBBv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
axb2.XLim = [bm, bM]; axb2.YLim = [-2, 6];
axb2.XTick = [bm, 0, 1, bM]; axb2.YTick = [];
axb2.Color = 'none';


pm = 1.32; pM = 1.52; p0 = 1.4; p1 = 1.44;
bm = (pm - p0) / (p1 - p0);
bM = (pM - p0) / (p1 - p0);

nexttile
plot($bfcnADh, $bfcnADv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
xlim([pm, pM]); ylim([-1, 15]); xticks([pm, p0, p1, pM]); yticks([-1, 15]);
grid on
text(0, 1.06, charlbl{3}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot($bfcnBDh, $bfcnBDv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
xlim([bm, bM]); ylim([-1, 15]); xticks([bm, 0, 1, bM]); yticks([-1, 15]);
grid on

axd1 = axes(fig);
axd1.Layout.Tile = 9;
plot($bfcnADh, $bfcnADv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
axd1.XLim = [pm, pM]; axd1.YLim = [-1, 15];
axd1.XTick = []; axd1.YTick = [-1, 15];

axd2 = axes(fig);
axd2.Layout.Tile = 9;
plot($bfcnBDh, $bfcnBDv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
axd2.XLim = [bm, bM]; axd2.YLim = [-1, 15];
axd2.XTick = [bm, 0, 1, bM]; axd2.YTick = [];
axd2.Color = 'none';
grid on

nexttile
plot($bfcnAFh, $bfcnAFv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
xlim([0.12, 0.16]); ylim([-0.5, 2]); xticks([0.12, 0.135, 0.14, 0.16]); yticks([-0.5, 2]);
grid on
text(0, 1.06, charlbl{4}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot($bfcnBFh, $bfcnBFv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
xlim([-3, 5]); ylim([-0.5, 2]); xticks([-3, 0, 1, 5]); yticks([-0.5, 2]);
grid on

axf1 = axes(fig);
axf1.Layout.Tile = 12;
plot($bfcnAFh, $bfcnAFv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
axf1.XLim = [0.12, 0.16]; axf1.YLim = [-0.5, 2];
axf1.XTick = []; axf1.YTick = [-0.5, 2];

axf2 = axes(fig);
axf2.Layout.Tile = 12;
plot($bfcnBFh, $bfcnBFv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
axf2.XLim = [-3, 5]; axf2.YLim = [-0.5, 2];
axf2.XTick = [-3, 0, 1, 5]; axf2.YTick = [-0.5, 2];
axf2.Color = 'none';
grid on


pm = 0.53; pM = 0.64; p0 = 0.575; p1 = 0.59;
bm = (pm - p0) / (p1 - p0);
bM = round((pM - p0) / (p1 - p0), 1);

nexttile
plot($bfcnARh, $bfcnARv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
xlim([pm, pM]); ylim([3, 6]); xticks([pm, p0, p1, pM]); yticks([3, 6]);
grid on
text(0, 1.06, charlbl{5}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot($bfcnBRh, $bfcnBRv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
xlim([bm, bM]); ylim([3, 6]); xticks([bm, 0, 1, bM]); yticks([3, 6]);
grid on

axr1 = axes(fig);
axr1.Layout.Tile = 15;
plot($bfcnARh, $bfcnARv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0])
axr1.XLim = [pm, pM]; axr1.YLim = [3, 6];
axr1.XTick = []; axr1.YTick = [3, 6];

axr2 = axes(fig);
axr2.Layout.Tile = 15;
plot($bfcnBRh, $bfcnBRv, 'o', 'MarkerSize', 0.2, 'LineStyle', 'none', 'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0])
axr2.XLim = [bm, bM]; axr2.YLim = [3, 6];
axr2.XTick = [bm, 0, 1, bM]; axr2.YTick = [3, 6];
axr2.Color = 'none';
grid on
"""

