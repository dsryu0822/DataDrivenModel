using MATLAB

trajA0_9596 = CSV.read("G:/BF/foodchain/trajA0_9596.csv", DataFrame)
trajA1_9596 = CSV.read("G:/BF/foodchain/trajA1_9596.csv", DataFrame)
trajB0_9596 = CSV.read("G:/BF/foodchain/trajB0_9596.csv", DataFrame)
trajB1_9596 = CSV.read("G:/BF/foodchain/trajB1_9596.csv", DataFrame)
trajC0_9596 = CSV.read("G:/BF/foodchain/trajC0_9596.csv", DataFrame)
trajC1_9596 = CSV.read("G:/BF/foodchain/trajC1_9596.csv", DataFrame)

bfcnAh, bfcnAv = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnA.jld2"))
bfcnBh, bfcnBv = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnB_9596.jld2"))
bfcnCh, bfcnCv = dict2bifurcation(callbfcn("G:/BF/foodchain/bfcnC_9596.jld2"))

mat"""
nIDs = 26;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

fig = tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
set(gcf, 'Position', [100, 100, 800, 600])

nexttile
plot3($(trajA0_9596.R), $(trajA0_9596.C), $(trajA0_9596.P), 'k', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = []);
text(0.01 ,0.95, charlbl{1}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3($(trajA1_9596.R), $(trajA1_9596.C), $(trajA1_9596.P), 'k', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
text(0.01 ,0.95, charlbl{2}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile([1 2])
scatter($bfcnAh, $bfcnAv, 1, '.k')
set(gca, xtick = [0.88, 0.95, 0.96, 1.0], ytick = [.55, .8])
grid on
box on
text(0.01 ,0.95, charlbl{7}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
text(0.27, -0.07, 'K', Units = 'normalized', FontSize = 12)
text(-0.04, .4,'P_{min}',Units = 'normalized', FontSize = 12, Rotation = 90)

nexttile
plot3($(trajB0_9596.R), $(trajB0_9596.C), $(trajB0_9596.P), 'r', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
text(0.01 ,0.95, charlbl{3}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3($(trajB1_9596.R), $(trajB1_9596.C), $(trajB1_9596.P), 'r', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
text(0.01 ,0.95, charlbl{4}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile([1 2])
scatter($bfcnBh, $bfcnBv, 1, '.r')
set(gca, xtick = [-7, 0, 1, 5], xlim = [-7, 5], ytick = [.55, .8])
grid on
box on
text(0.01 ,0.95, charlbl{8}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
text(0.27, -0.07, 'K', Units = 'normalized', FontSize = 12)
text(-0.04, .4,'P_{min}',Units = 'normalized', FontSize = 12, Rotation = 90)

nexttile
plot3($(trajC0_9596.R), $(trajC0_9596.C), $(trajC0_9596.P), 'b', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
text(0.01 ,0.95, charlbl{5}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3($(trajC1_9596.R), $(trajC1_9596.C), $(trajC1_9596.P), 'b', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
text(0.01 ,0.95, charlbl{6}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile([1 2])
scatter($bfcnCh, $bfcnCv, 1, '.b')
set(gca, xtick = [-7, 0, 1, 5], xlim = [-7, 5], ytick = [.55, .8])
grid on
box on
text(0.01 ,0.95, charlbl{9}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
text(0.27, -0.07, 'K', Units = 'normalized', FontSize = 12)
text(-0.04, .4,'P_{min}',Units = 'normalized', FontSize = 12, Rotation = 90)
"""
