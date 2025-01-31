close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               soft impact                     %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = readtable("G:/DDM/result/soft_shortterm.csv");
subplot(3, 1, 1)
hold on
plot(dat.t_t, dat.t_u, '-', 'LineWidth', 2, 'Color', [0, 0, 0, .5])
plot(dat.p_t, dat.p_u, '-', 'LineWidth', 2, 'Color', [0, 0, 1, .5])
ylabel('$u$', 'interpreter', 'latex', 'FontSize', 13);
xlim([0, 12]); xticks([0, 12]); yticks([])
% legend('ground truth', 'recovered', 'FontSize', 14)
xlabel0 = xlabel('$t$', 'interpreter', 'latex', 'FontSize', 12);
xlabel0.Position(2) = xlabel0.Position(2) + 0.01;
hold off
text(0.0,1.1, '(a)','Units','normalized');
set(gca, 'Position', [0.07 0.7 0.9 0.23])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               gear system                     %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = readtable("G:/DDM/result/gear_shortterm.csv");
subplot(3, 1, 2)
hold on
plot(dat.t__, dat.t_x, '-', 'LineWidth', 2, 'Color', [0, 0, 0, .5])
plot(dat.p__, dat.p_x, '-', 'LineWidth', 2, 'Color', [0, 0, 1, .5])
ylabel('$x$', 'interpreter', 'latex', 'FontSize', 13);
xticks([0, 250]); yticks([]); xlim([0, 250]);
xlabel0 = xlabel('$t$', 'interpreter', 'latex', 'FontSize', 12);
xlabel0.Position(2) = xlabel0.Position(2) + .4;
hold off
text(0.0,1.1, '(b)','Units','normalized');
set(gca, 'Position', [0.07 0.38 0.9 0.23])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               Hindmarsh-Rose                  %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = readtable("G:/DDM/result/hrnm_shortterm.csv");
subplot(3, 1, 3)
hold on
plot(dat.t_t, dat.t_z, '-', 'LineWidth', 2, 'Color', [0, 0, 0, .5])
plot(dat.p_t, dat.p_z, '-', 'LineWidth', 2, 'Color', [0, 0, 1, .5])
ylabel('$z$', 'interpreter', 'latex', 'FontSize', 13);
xticks([0, 200]); yticks([]);
xlabel0 = xlabel('$t$', 'interpreter', 'latex', 'FontSize', 12);
xlabel0.Position(2) = xlabel0.Position(2) + 0.9;
hold off
text(0.0,1.1, '(c)','Units','normalized');
set(gca, 'Position', [0.07 0.06 0.9 0.23])

set(gcf,'Position', [3000 100 400 400])
