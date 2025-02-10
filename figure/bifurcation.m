%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               soft impact                     %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = readtable("G:/DDM/result/soft_bfcn.csv");
subplot(2, 3, 1)
scatter(dat.hrzn, dat.vrtc, 1, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'k')
xticks([]);
yticks([]); ylabel('$\dot{u}$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'Position', [0.03 0.51 0.3 0.45])
text(0.0,1.05, '(a)','Units','normalized', 'FontSize', 14)

dat = readtable("G:/DDM/result/soft_bfcn_rcvd.csv");
subplot(2, 3, 4)
scatter(dat.hrzn, dat.vrtc, 1, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'b')
xticks([min(dat.hrzn), max(dat.hrzn)]);
yticks([]); ylabel('$\dot{u}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0 = xlabel('$d$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0.Position(2) = xlabel0.Position(2) + 0.08;
set(gca, 'Position', [0.03 0.05 0.3 0.45])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               gear system                     %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = readtable("G:/DDM/result/gear_bfcn.csv");
subplot(2, 3, 2)
scatter(dat.hrzn, dat.vrtc, 1, 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerFaceColor', 'k')
xticks([]);
yticks([]); ylabel('$v$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'Position', [0.36 0.51 0.3 0.45])
text(0.0,1.05, '(b)','Units','normalized', 'FontSize', 14)

dat = readtable("G:/DDM/result/gear_bfcn_rcvd.csv");
subplot(2, 3, 5)
scatter(dat.hrzn, dat.vrtc, 1, 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerFaceColor', 'b')
xticks([min(dat.hrzn), max(dat.hrzn)]);
yticks([]); ylabel('$v$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0 = xlabel('$\bar{F}_{e}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0.Position(2) = xlabel0.Position(2) + 0.35;
set(gca, 'Position', [0.36 0.05 0.3 0.45])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               Hindmarsh-Rose                  %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = readtable("G:/DDM/result/hrnm_bfcn.csv");
subplot(2, 3, 3)
scatter(dat.hrzn, dat.vrtc, 1, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'k')
xticks([]);
yticks([]); ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'Position', [0.69 0.51 0.3 0.45])
text(0.0,1.05, '(c)','Units','normalized', 'FontSize', 14)

dat = readtable("G:/DDM/result/hrnm_bfcn_rcvd.csv");
subplot(2, 3, 6)
scatter(dat.hrzn, dat.vrtc, 1, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'b')
xticks([min(dat.hrzn), max(dat.hrzn)]);
yticks([]); ylabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0 = xlabel('$l$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0.Position(2) = xlabel0.Position(2) + 0.4;
set(gca, 'Position', [0.69 0.05 0.3 0.45])

set(gcf,'Position', [3000 100 1000 600])
