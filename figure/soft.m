close all

shot_term = readtable("G:/DDM/result/soft_shortterm.csv");
subplot(5, 1, 1)
plot(shot_term.t_t, shot_term.t_u, 'LineWidth', 1.5, 'Color', [0, 0, 0, .5]);
hold on
plot(shot_term.t_t, shot_term.p_u, 'LineWidth', 1.5, 'Color', [0, 0, 1, .5]);
xlim([0, 12]); xticks([0, 12]); xlabel('$t$', 'interpreter', 'latex', 'FontSize', 12);
yticks([]); ylabel('$u$', 'interpreter', 'latex', 'FontSize', 13);
text(0.0,1.06, '(a)','Units','normalized', 'FontSize', 14)
hold off

bifurcation = readtable("G:/DDM/result/soft_bfcn.csv");
subplot(5, 1, 2)
scatter(bifurcation.hrzn, bifurcation.vrtc, 1, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor', 'k')
xticks([0.1, 0.3]); xlabel('$d$', 'interpreter', 'latex', 'FontSize', 12);
yticks([]); ylabel('$v$', 'Interpreter', 'latex', 'FontSize', 14);
text(0.0,1.06, '(b)','Units','normalized', 'FontSize', 14)
bifurcation = readtable("G:/DDM/result/soft_bfcn_rcvd.csv");
subplot(5, 1, 3)
scatter(bifurcation.hrzn, bifurcation.vrtc, 1, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerFaceColor', 'b')
xticks([0.1, 0.3]); xlabel('$d$', 'interpreter', 'latex', 'FontSize', 12);
yticks([]); ylabel('$v$', 'Interpreter', 'latex', 'FontSize', 14);
text(0.0,1.06, '(c)','Units','normalized', 'FontSize', 14)

lyapunov = readtable("G:/DDM/result/soft_lyapunov.csv");
subplot(5, 1, 4)
hold on
plot(lyapunov.bp, lyapunov.x_1, 'LineWidth', 1.5, 'Color', [0, 0, 0, .5])
plot(lyapunov.bp, lyapunov.x_2, 'LineWidth', 1.5, 'Color', [0, 0, 0, .5])
plot(lyapunov.bp, lyapunov.x_3, 'LineWidth', 1.5, 'Color', [0, 0, 0, .5])
xticks([0.1, 0.3]); xlabel('$d$', 'interpreter', 'latex', 'FontSize', 12);
yticks([-3, 0, 1]); ylabel('$\lambda$', 'interpreter', 'latex', 'FontSize', 13); 
text(0.0,1.06, '(d)','Units','normalized', 'FontSize', 14)
hold off
lyapunov = readtable("G:/DDM/result/soft_lyapunov_rcvd.csv");
subplot(5, 1, 5)
hold on
plot(lyapunov.bp, lyapunov.x_1, 'LineWidth', 1.5, 'Color', [0, 0, 1, .5])
plot(lyapunov.bp, lyapunov.x_2, 'LineWidth', 1.5, 'Color', [0, 0, 1, .5])
plot(lyapunov.bp, lyapunov.x_3, 'LineWidth', 1.5, 'Color', [0, 0, 1, .5])
xticks([0.1, 0.3]); xlabel('$d$', 'interpreter', 'latex', 'FontSize', 12);
yticks([-3, 0, 1]); ylabel('$\lambda$', 'interpreter', 'latex', 'FontSize', 13); 
text(0.0,1.06, '(e)','Units','normalized', 'FontSize', 14)
hold off

set(gcf,'Position', [3000 100 800 1600])