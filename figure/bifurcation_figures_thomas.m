set(0,'DefaultAxesFontname', 'Arial')
set(0,'DefaultTextFontname', 'Arial')


nIDs = 26;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

traj1 = readtable("traj1.csv"); head(traj1)
traj2 = readtable("traj2.csv"); head(traj2)
fail1 = readtable("fail1.csv"); head(fail1)
fail2 = readtable("fail2.csv"); head(fail2)
rcvd1 = readtable("rcvd1.csv"); head(rcvd1)
rcvd2 = readtable("rcvd2.csv"); head(rcvd2)
bfcn = readtable("bfcn.csv"); head(bfcn)
bfcn_fail = readtable("bfcn__.csv"); head(bfcn_fail)
bfcn_rcvd = readtable("bfcn_.csv"); head(bfcn_rcvd)

%%
clf
plt = tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
set(gcf, 'Position', [100, 100, 800, 600])

nexttile
plot3(traj1.x, traj1.y, traj1.z, 'k', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
zlabel('ground truth', FontSize = 12)
box on
text(0.01,0.93,charlbl{1}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3(traj2.x, traj2.y, traj2.z, 'k', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
box on
text(0.01,0.93,charlbl{2}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile([1, 2])
scatter(bfcn.h, bfcn.v, 1, '.k')
set(gca, xtick = [0.14, 0.141, 0.16], ytick = [-5, 7], ylim = [-5, 7])
grid on
box on
xtickangle(45)
text(0.01,0.93,charlbl{7}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
text(0.5, -0.08,'b',Units = 'normalized', FontSize = 12)
text(-0.04, .4,'x_{max}',Units = 'normalized', FontSize = 12, Rotation = 90)

nexttile
plot3(fail1.x1, fail1.x2, fail1.x3, 'r', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
zlabel('failed', FontSize = 12)
box on
text(0.01,0.93,charlbl{3}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3(fail2.x1, fail2.x2, fail2.x3, 'r', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
box on
text(0.01,0.93,charlbl{4}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile([1, 2])
scatter(bfcn_fail.h, bfcn_fail.v, 1, '.r')
set(gca, xtick = [0, 1, 20], ytick = [-5, 7], ylim = [-5, 7])
grid on
box on
text(0.01,0.93,charlbl{8}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
text(0.5, -0.08,'\beta',Units = 'normalized', FontSize = 12)
text(-0.04, .4,'x_{max}',Units = 'normalized', FontSize = 12, Rotation = 90)

nexttile
plot3(rcvd1.x1, rcvd1.x2, rcvd1.x3, 'b', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
zlabel('recovered', FontSize = 12)
box on
text(0.01,0.93,charlbl{5}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile
plot3(rcvd2.x1, rcvd2.x2, rcvd2.x3, 'b', 'LineWidth', 1)
set(gca, xtick = [], ytick = [], ztick = [])
box on
text(0.01,0.93,charlbl{6}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')

nexttile([1, 2])
scatter(bfcn_rcvd.h, bfcn_rcvd.v, 1, '.b')
set(gca, xtick = [0, 1, 20], ytick = [-5, 7], ylim = [-5, 7])
grid on
box on
text(0.01,0.93,charlbl{9}, Units = 'normalized', FontSize = 12, FontWeight = 'bold')
text(0.5, -0.08,'\beta',Units = 'normalized', FontSize = 12)
text(-0.04, .4,'x_{max}',Units = 'normalized', FontSize = 12, Rotation = 90)
savefig("thomas.fig")