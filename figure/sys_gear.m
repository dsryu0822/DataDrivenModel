close all

x = .1;
y = .195;
w = .85;
h = .13;

shot_term = readtable("G:/DDM/result/gear_shortterm.csv");
subplot(5, 1, 1)
set(gca, 'Position', [x 0.05+4*y w h]);
plot(shot_term.t__, shot_term.t_x, 'Color', [0, 0, 1, .5]);
hold on
plot(shot_term.t__, shot_term.p_x, 'Color', [1, 0, 0, .5]);
xlim([0, 250]); xticks([0, 125, 250]);
ylim([-0.5, 2.5]); yticks([0, 1, 2]);
xlabel('$t$', 'interpreter', 'latex', 'Position', [155, -0.8]);
ylabel('$x$', 'interpreter', 'latex', 'Position', [-10, 1], 'FontSize', 12);
text(0.0, 1.15, '(a)','Units','normalized')
hold off
set(gca, 'TickLength', [0, 0]);

bifurcation = readtable("G:/DDM/result/gear_bfcn.csv");
bifurcation = bifurcation(1:5:end, :);
subplot(5, 1, 2)
set(gca, 'Position', [x 0.05+3*y w h])
scatter(bifurcation.hrzn, bifurcation.vrtc, .5, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'k')
box on
xticks([-.2, 0, .2]); yticks([-2, 0, 2]);
ylim([-2.1, 2.1]);
xlabel('$\bar{F}_{e}$', 'interpreter', 'latex', 'Position', [0.05, -2.5]);
ylabel('$v$', 'Interpreter', 'latex', 'Position', [-.22, 0], 'FontSize', 12);
text(0.0, 1.15, '(b)','Units','normalized')
bifurcation = readtable("G:/DDM/result/gear_bfcn_rcvd.csv");
set(gca, 'TickLength', [0, 0]);
subplot(5, 1, 3)
set(gca, 'Position', [x 0.05+2*y w h])
scatter(bifurcation.hrzn, bifurcation.vrtc, .5, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'b')
box on
xticks([-.2, 0, .2]); yticks([-2, 0, 2]);
ylim([-2.1, 2.1]);
xlabel('$\bar{F}_{e}$', 'interpreter', 'latex', 'Position', [0.05, -2.5]);
ylabel('$v$', 'Interpreter', 'latex', 'Position', [-.22, 0], 'FontSize', 12);
text(0.0, 1.15, '(c)','Units','normalized')
set(gca, 'TickLength', [0, 0]);

lyapunov = readtable("G:/DDM/result/gear_lyapunov.csv");
subplot(5, 1, 4)
set(gca, 'Position', [x 0.05+1*y w h])
hold on
plot(lyapunov.bp, lyapunov.x_1, 'Color', [0, 0, 0, .5])
plot(lyapunov.bp, lyapunov.x_2, 'Color', [0, 0, 0, .5])
plot(lyapunov.bp, lyapunov.x_3, 'Color', [0, 0, 0, .5])
box on
xticks([-.2, 0, .2]);
ylim([-.07, 0.04]); yticks([-.06, 0, .03]);
xlabel('$\bar{F}_{e}$', 'interpreter', 'latex', 'Position', [0.05, -.08]);
ylabel('$\lambda$', 'interpreter', 'latex', 'Position', [-.22, -0.02]);
text(0.0, 1.15, '(d)','Units','normalized')
set(gca, 'TickLength', [0, 0]);
hold off
lyapunov = readtable("G:/DDM/result/gear_lyapunov_rcvd.csv");
subplot(5, 1, 5)
set(gca, 'Position', [x 0.05+0*y w h])
hold on
plot(lyapunov.bp, lyapunov.x_1, 'Color', [0, 0, 1, .5])
plot(lyapunov.bp, lyapunov.x_2, 'Color', [0, 0, 1, .5])
plot(lyapunov.bp, lyapunov.x_3, 'Color', [0, 0, 1, .5])
box on
xticks([-.2, 0, .2]);
ylim([-.07, 0.04]); yticks([-.06, 0, .03]);
xlabel('$\bar{F}_{e}$', 'interpreter', 'latex', 'Position', [0.05, -.08]);
ylabel('$\lambda$', 'interpreter', 'latex', 'Position', [-.22, -0.02]);
text(0.0, 1.15, '(e)','Units','normalized')
hold off
set(gca, 'TickLength', [0, 0]);

set(gcf, 'Position', [3100 100 300 600])