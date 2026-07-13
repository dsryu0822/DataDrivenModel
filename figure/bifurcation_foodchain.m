nIDs = 26;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

trajA0 = readtable("G:/BF/foodchain/trajA0.csv");
trajA1 = readtable("G:/BF/foodchain/trajA1.csv");
trajB0 = readtable("G:/BF/foodchain/trajB0.csv");
trajB1 = readtable("G:/BF/foodchain/trajB1.csv");
trajC0 = readtable("G:/BF/foodchain/trajC0.csv");
trajC1 = readtable("G:/BF/foodchain/trajC1.csv");
trajA0 = trajA0(1:100:end, :);
trajA1 = trajA1(1:100:end, :);

%%
clf
fig = tiledlayout(2,2);
set(gcf, 'Position', [100, 100, 500, 500])

nexttile
plot3(trajA0.R, trajA0.C, trajA0.P, 'o', 'MarkerSize', 3, 'Color', [0.5, 0.5, 0.5])
hold on
plot3(trajB0.R, trajB0.C, trajB0.P, 'r', 'LineWidth', 1)
text(0.01,0.93,charlbl{1}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3(trajA1.R, trajA1.C, trajA1.P, 'o', 'MarkerSize', 3, 'Color', [0.5, 0.5, 0.5])
hold on
plot3(trajB1.R, trajB1.C, trajB1.P, 'r', 'LineWidth', 1)
text(0.01,0.93,charlbl{2}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3(trajA0.R, trajA0.C, trajA0.P, 'o', 'MarkerSize', 3, 'Color', [0.5, 0.5, 0.5])
hold on
plot3(trajC0.R, trajC0.C, trajC0.P, 'b', 'LineWidth', 1)
text(0.01,0.93,charlbl{3}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3(trajA1.R, trajA1.C, trajA1.P, 'o', 'MarkerSize', 3, 'Color', [0.5, 0.5, 0.5])
hold on
plot3(trajC1.R, trajC1.C, trajC1.P, 'b', 'LineWidth', 1)
text(0.01,0.93,charlbl{4}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
