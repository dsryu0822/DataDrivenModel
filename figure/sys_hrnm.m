close all

x = .1;
y = .195;
w = .85;
h = .13;

shot_term = readtable("G:/DDM/result/hrnm_shortterm.csv");
subplot(5, 1, 1)
set(gca, 'Position', [x 0.05+4*y w h]);
plot(shot_term.t_t, shot_term.t_x, 'Color', [0, 0, 1, .5]);
hold on
plot(shot_term.t_t, shot_term.p_x, 'Color', [1, 0, 0, .5]);
xlim([0, 200]); xticks([0, 100, 200]);
ylim([-3, 3]); yticks([-2, 0, 2]);
xlabel('$t$', 'interpreter', 'latex', 'Position', [120, -3.5]);
ylabel('$x$', 'interpreter', 'latex', 'Position', [-10, 0], 'FontSize', 12);
text(0.0, 1.15, '(a)','Units','normalized')
hold off
set(gca, 'TickLength', [0, 0]);

bifurcation = readtable("G:/DDM/result/hrnm_bfcn.csv");
subplot(5, 1, 2)
set(gca, 'Position', [x 0.05+3*y w h])
scatter(bifurcation.hrzn, bifurcation.vrtc, .5, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'k')
box on
xticks([0, 0.5, 1]); yticks([-2, 0, 2]);
ylim([-2.3, 2.5]);
xlabel('$l$', 'interpreter', 'latex', 'Position', [0.6, -2.7]);
ylabel('$x$', 'Interpreter', 'latex', 'Position', [-.05, 0], 'FontSize', 12);
text(0.0, 1.15, '(b)','Units','normalized')
bifurcation = readtable("G:/DDM/result/hrnm_bfcn_rcvd.csv");
set(gca, 'TickLength', [0, 0]);
subplot(5, 1, 3)
set(gca, 'Position', [x 0.05+2*y w h])
scatter(bifurcation.hrzn, bifurcation.vrtc, .5, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'b')
box on
xticks([0, 0.5, 1]); yticks([-2, 0, 2]);
ylim([-2.3, 2.5]);
xlabel('$l$', 'interpreter', 'latex', 'Position', [0.6, -2.7]);
ylabel('$x$', 'Interpreter', 'latex', 'Position', [-.05, 0], 'FontSize', 12);
text(0.0, 1.15, '(c)','Units','normalized')
set(gca, 'TickLength', [0, 0]);

lyapunov = readtable("G:/DDM/result/hrnm_lyapunov.csv");
subplot(5, 1, 4)
set(gca, 'Position', [x 0.05+1*y w h])
hold on
plot(lyapunov.bp, lyapunov.x_1, 'Color', [0, 0, 0, .5])
plot(lyapunov.bp, lyapunov.x_2, 'Color', [0, 0, 0, .5])
plot(lyapunov.bp, lyapunov.x_3, 'Color', [0, 0, 0, .5])
box on
xticks([0, 0.5, 1]);
ylim([-.7, .2]); yticks([-.6, 0, .2]);
xlabel('$l$', 'interpreter', 'latex', 'Position', [0.6, -0.8]);
ylabel('$\lambda$', 'interpreter', 'latex', 'Position', [-.05, -.3]);
text(0.0, 1.15, '(d)','Units','normalized')
set(gca, 'TickLength', [0, 0]);
hold off
lyapunov = readtable("G:/DDM/result/hrnm_lyapunov_rcvd.csv");
subplot(5, 1, 5)
set(gca, 'Position', [x 0.05+0*y w h])
hold on
plot(lyapunov.bp, lyapunov.x_1, 'Color', [0, 0, 1, .5])
plot(lyapunov.bp, lyapunov.x_2, 'Color', [0, 0, 1, .5])
plot(lyapunov.bp, lyapunov.x_3, 'Color', [0, 0, 1, .5])
box on
xticks([0, 0.5, 1]);
ylim([-.7, .2]); yticks([-.6, 0, .2]);
xlabel('$l$', 'interpreter', 'latex', 'Position', [0.6, -0.8]);
ylabel('$\lambda$', 'interpreter', 'latex', 'Position', [-.05, -.3]);
text(0.0, 1.15, '(e)','Units','normalized')
hold off
set(gca, 'TickLength', [0, 0]);

set(gcf, 'Position', [3100 100 300 600])