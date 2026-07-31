using MATLAB

bfcnAh, bfcnAv = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnA.jld2"));

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
scatter($bfcnAh, $bfcnAv, 0.5, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('0.93   0.94    ', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axa1.XAxisLocation = 'top';
axa1.XLim = [0.88, 1.00]; axa1.YLim = [0.55, 0.8];
axa1.XTick = []; axa1.YTick = [0.55, 0.8];
axa2 = axes(fig);
axa2.Layout.Tile = 1;
scatter($bfcnB_9394h, $bfcnB_9394v, 1, '.r')
box on
grid on
axa2.XLim = [-5, 7]; axa2.YLim = [0.55, 0.8];
axa2.XTick = [-5, 0, 1, 7]; axa2.YTick = [];
axa2.Color = 'none';
text(0.02 , 1.05, charlbl{1}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axb1 = axes(fig);
axb1.Layout.Tile = 2;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('0.93     0.95', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axb1.XAxisLocation = 'top';
axb1.XLim = [0.88, 1.00]; axb1.YLim = [0.55, 0.8];
axb1.XTick = []; axb1.YTick = [0.55, 0.8];
axb2 = axes(fig);
axb2.Layout.Tile = 2;
scatter($bfcnB_9395h, $bfcnB_9395v, 1, '.r')
box on
grid on
axb2.XLim = [-2.5, 3.5]; axb2.YLim = [0.55, 0.8];
axb2.XTick = [-2.5, 0, 1, 3.5]; axb2.YTick = [];
axb2.Color = 'none';
text(0.02 , 1.05, charlbl{2}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axc1 = axes(fig);
axc1.Layout.Tile = 3;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('   0.94   0.95', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axc1.XAxisLocation = 'top';
axc1.XLim = [0.88, 1.00]; axc1.YLim = [0.55, 0.8];
axc1.XTick = []; axc1.YTick = [0.55, 0.8];
axc2 = axes(fig);
axc2.Layout.Tile = 3;
scatter($bfcnB_9495h, $bfcnB_9495v, 1, '.r')
box on
grid on
axc2.XLim = [-6, 6]; axc2.YLim = [0.55, 0.8];
axc2.XTick = [-6, 0, 1, 6]; axc2.YTick = [];
axc2.Color = 'none';
text(0.02 , 1.05, charlbl{3}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axd1 = axes(fig);
axd1.Layout.Tile = 4;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('   0.93        0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axd1.XAxisLocation = 'top';
axd1.XLim = [0.88, 1.00]; axd1.YLim = [0.55, 0.8];
axd1.XTick = []; axd1.YTick = [0.55, 0.8];
axd2 = axes(fig);
axd2.Layout.Tile = 4;
scatter($bfcnB_9396h, $bfcnB_9396v, 1, '.r')
box on
grid on
axd2.XLim = [-1.7, 2.4]; axd2.YLim = [0.55, 0.8];
axd2.XTick = [-1.7, 0, 1, 2.4]; axd2.YTick = [];
axd2.Color = 'none';
text(0.02 , 1.05, charlbl{4}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axe1 = axes(fig);
axe1.Layout.Tile = 5;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('         0.94   0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axe1.XAxisLocation = 'top';
axe1.XLim = [0.88, 1.00]; axe1.YLim = [0.55, 0.8];
axe1.XTick = []; axe1.YTick = [0.55, 0.8];
axe2 = axes(fig);
axe2.Layout.Tile = 5;
scatter($bfcnB_9496h, $bfcnB_9496v, 1, '.r')
box on
grid on
axe2.XLim = [-3, 3]; axe2.YLim = [0.55, 0.8];
axe2.XTick = [-3, 0, 1, 3]; axe2.YTick = [];
axe2.Color = 'none';
text(0.02 , 1.05, charlbl{5}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axf1 = axes(fig);
axf1.Layout.Tile = 6;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('               0.95   0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axf1.XAxisLocation = 'top';
axf1.XLim = [0.88, 1.00]; axf1.YLim = [0.55, 0.8];
axf1.XTick = []; axf1.YTick = [0.55, 0.8];
axf2 = axes(fig);
axf2.Layout.Tile = 6;
scatter($bfcnB_9596h, $bfcnB_9596v, 1, '.r')
box on
grid on
axf2.XLim = [-7, 5]; axf2.YLim = [0.55, 0.8];
axf2.XTick = [-7, 0, 1, 5]; axf2.YTick = [];
axf2.Color = 'none';
text(0.02 , 1.05, charlbl{6}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axg1 = axes(fig);
axg1.Layout.Tile = 7;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('0.93   0.94    ', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axg1.XAxisLocation = 'top';
axg1.XLim = [0.88, 1.00]; axg1.YLim = [0.55, 0.8];
axg1.XTick = []; axg1.YTick = [0.55, 0.8];
axg2 = axes(fig);
axg2.Layout.Tile = 7;
scatter($bfcnC_9394h, $bfcnC_9394v, 1, '.b')
box on
grid on
axg2.XLim = [-5, 7]; axg2.YLim = [0.55, 0.8];
axg2.XTick = [-5, 0, 1, 7]; axg2.YTick = [];
axg2.Color = 'none';
text(0.02 , 1.05, charlbl{7}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axh1 = axes(fig);
axh1.Layout.Tile = 8;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('0.93     0.95', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axh1.XAxisLocation = 'top';
axh1.XLim = [0.88, 1.00]; axh1.YLim = [0.55, 0.8];
axh1.XTick = []; axh1.YTick = [0.55, 0.8];
axh2 = axes(fig);
axh2.Layout.Tile = 8;
scatter($bfcnC_9395h, $bfcnC_9395v, 1, '.b')
box on
grid on
axh2.XLim = [-2.5, 3.5]; axh2.YLim = [0.55, 0.8];
axh2.XTick = [-2.5, 0, 1, 3.5]; axh2.YTick = [];
axh2.Color = 'none';
text(0.02 , 1.05, charlbl{8}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axi1 = axes(fig);
axi1.Layout.Tile = 9;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('   0.94   0.95', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axi1.XAxisLocation = 'top';
axi1.XLim = [0.88, 1.00]; axi1.YLim = [0.55, 0.8];
axi1.XTick = []; axi1.YTick = [0.55, 0.8];
axi2 = axes(fig);
axi2.Layout.Tile = 9;
scatter($bfcnC_9495h, $bfcnC_9495v, 1, '.b')
box on
grid on
axi2.XLim = [-6, 6]; axi2.YLim = [0.55, 0.8];
axi2.XTick = [-6, 0, 1, 6]; axi2.YTick = [];
axi2.Color = 'none';
text(0.02 , 1.05, charlbl{9}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axj1 = axes(fig);
axj1.Layout.Tile = 10;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('   0.93        0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axj1.XAxisLocation = 'top';
axj1.XLim = [0.88, 1.00]; axj1.YLim = [0.55, 0.8];
axj1.XTick = []; axj1.YTick = [0.55, 0.8];
axj2 = axes(fig);
axj2.Layout.Tile = 10;
scatter($bfcnC_9396h, $bfcnC_9396v, 1, '.b')
box on
grid on
axj2.XLim = [-1.7, 2.4]; axj2.YLim = [0.55, 0.8];
axj2.XTick = [-1.7, 0, 1, 2.4]; axj2.YTick = [];
axj2.Color = 'none';
text(0.02 , 1.05, charlbl{10}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axk1 = axes(fig);
axk1.Layout.Tile = 11;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('         0.94   0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axk1.XAxisLocation = 'top';
axk1.XLim = [0.88, 1.00]; axk1.YLim = [0.55, 0.8];
axk1.XTick = []; axk1.YTick = [0.55, 0.8];
axk2 = axes(fig);
axk2.Layout.Tile = 11;
scatter($bfcnC_9496h, $bfcnC_9496v, 1, '.b')
box on
grid on
axk2.XLim = [-3, 3]; axk2.YLim = [0.55, 0.8];
axk2.XTick = [-3, 0, 1, 3]; axk2.YTick = [];
axk2.Color = 'none';
text(0.02 , 1.05, charlbl{11}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

axl1 = axes(fig);
axl1.Layout.Tile = 12;
scatter($bfcnAh, $bfcnAv, 1, '.', 'MarkerEdgeColor', [0.7, 0.7, 0.7])
title('               0.95   0.96', FontSize = 10)
text(0.2, -0.1, '\\beta', Units = 'normalized', FontSize = 12)
text(-0.07, 0.5, 'P', Units = 'normalized', FontSize = 12, Rotation = 90)
axl1.XAxisLocation = 'top';
axl1.XLim = [0.88, 1.00]; axl1.YLim = [0.55, 0.8];
axl1.XTick = []; axl1.YTick = [0.55, 0.8];
axl2 = axes(fig);
axl2.Layout.Tile = 12;
scatter($bfcnC_9596h, $bfcnC_9596v, 1, '.b')
box on
grid on
axl2.XLim = [-7, 5]; axl2.YLim = [0.55, 0.8];
axl2.XTick = [-7, 0, 1, 5]; axl2.YTick = [];
axl2.Color = 'none';
text(0.02 , 1.05, charlbl{12}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
"""


# nexttile
# scatter($bfcnC_9396h, $bfcnC_9396v, 1, '.b')

# nexttile
# scatter($bfcnC_9495h, $bfcnC_9495v, 1, '.b')

# nexttile
# scatter($bfcnC_9496h, $bfcnC_9496v, 1, '.b')

# nexttile
# scatter($bfcnC_9596h, $bfcnC_9596v, 1, '.b')
