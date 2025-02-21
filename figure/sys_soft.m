close all

x = .1;
y = .195;
w = .85;
h = .13;

shot_term = readtable("G:/DDM/result/soft_shortterm.csv");
subplot(5, 1, 1)
set(gca, 'Position', [x 0.05+4*y w h]);
plot(shot_term.t_t, shot_term.t_u, 'Color', [0, 0, 1, .5]);
hold on
plot(shot_term.t_t, shot_term.p_u, 'Color', [1, 0, 0, .5]);
xlim([0, 12]); xticks([0, 6, 12]);
ylim([-0.07, 0.07]); yticks([-0.05, 0, 0.05]);
xlabel('$t$', 'interpreter', 'latex', 'Position', [7, -0.085]);
ylabel('$u$', 'interpreter', 'latex', 'Position', [-.5, 0], 'FontSize', 12);
text(0.0, 1.15, '(a)','Units','normalized')
hold off
set(gca, 'TickLength', [0, 0]);

bifurcation = readtable("G:/DDM/result/soft_bfcn.csv");
subplot(5, 1, 2)
set(gca, 'Position', [x 0.05+3*y w h])
scatter(bifurcation.hrzn, bifurcation.vrtc, .5, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'k')
box on
xticks([0.1, 0.2, 0.3]); yticks([-0.3, 0, 0.3]);
ylim([-.4, .4]);
xlabel('$d$', 'interpreter', 'latex', 'Position', [0.22, -0.5]);
ylabel('$v$', 'Interpreter', 'latex', 'Position', [0.092, 0], 'FontSize', 12);
text(0.0, 1.15, '(b)','Units','normalized')
bifurcation = readtable("G:/DDM/result/soft_bfcn_rcvd.csv");
bifurcation = bifurcation(1:2:end, :);
set(gca, 'TickLength', [0, 0]);
subplot(5, 1, 3)
set(gca, 'Position', [x 0.05+2*y w h])
scatter(bifurcation.hrzn, bifurcation.vrtc, .5, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'b')
box on
xticks([0.1, 0.2, 0.3]); yticks([-0.3, 0, 0.3]);
ylim([-.4, .4]);
xlabel('$d$', 'interpreter', 'latex', 'Position', [0.22, -0.5]);
ylabel('$v$', 'Interpreter', 'latex', 'Position', [0.092, 0], 'FontSize', 12);
text(0.0, 1.15, '(c)','Units','normalized')
set(gca, 'TickLength', [0, 0]);

lyapunov = readtable("G:/DDM/result/soft_lyapunov.csv");
subplot(5, 1, 4)
set(gca, 'Position', [x 0.05+1*y w h])
hold on
plot(lyapunov.bp, lyapunov.x_1, 'Color', [0, 0, 0, .5])
plot(lyapunov.bp, lyapunov.x_2, 'Color', [0, 0, 0, .5])
plot(lyapunov.bp, lyapunov.x_3, 'Color', [0, 0, 0, .5])
box on
xticks([0.1, 0.2, 0.3]);
ylim([-3.2, 1.2]); yticks([-3, 0, 1]);
xlabel('$d$', 'interpreter', 'latex', 'Position', [0.22, -3.7]);
ylabel('$\lambda$', 'interpreter', 'latex', 'Position', [0.092, -1]);
text(0.0, 1.15, '(d)','Units','normalized')
set(gca, 'TickLength', [0, 0]);
hold off
lyapunov = readtable("G:/DDM/result/soft_lyapunov_rcvd.csv");
subplot(5, 1, 5)
set(gca, 'Position', [x 0.05+0*y w h])
hold on
plot(lyapunov.bp, lyapunov.x_1, 'Color', [0, 0, 1, .5])
plot(lyapunov.bp, lyapunov.x_2, 'Color', [0, 0, 1, .5])
plot(lyapunov.bp, lyapunov.x_3, 'Color', [0, 0, 1, .5])
box on
xticks([0.1, 0.2, 0.3]);
ylim([-3.2, 1.2]); yticks([-3, 0, 1]);
xlabel('$d$', 'interpreter', 'latex', 'Position', [0.22, -3.7]);
ylabel('$\lambda$', 'interpreter', 'latex', 'Position', [0.092, -1]);
text(0.0, 1.15, '(e)','Units','normalized')
hold off
set(gca, 'TickLength', [0, 0]);

set(gcf, 'Position', [3100 100 300 600])