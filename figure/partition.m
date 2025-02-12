dat = readtable("G:/DDM/_figure/partition.csv");

x = 0.1;
w = 0.85;
h = 0.24;
subplot(3,1,1)
set(gca, 'Position', [x 0.71 w h])
p1 = plot(dat.softh, 100*dat.softp, 'k.');
ylabel('detected (%)', 'Position', [8e-6, 50])
xscale log; xlabel('$h$','interpreter','latex', 'Position', [3e-4, -20]);
xlim([.9999e-5, 1e-2]); ylim([-10, 110]); yticks([0, 100]);
text(0.0,1.1, '(a)','Units','normalized', 'FontSize', 12)

subplot(3,1,2)
set(gca, 'Position', [x 0.38 w h])
p2 = plot(dat.gearh, 100*dat.gearp, 'k.');
ylabel('detected (%)', 'Position', [8e-4, 50])
xscale log; xlabel('$h$','interpreter','latex', 'Position', [3e-2, -20]);
xlim([1e-3, 1e-0]); ylim([-10, 110]); yticks([0, 100]);
text(0.0,1.1, '(b)','Units','normalized', 'FontSize', 12)

subplot(3,1,3)
set(gca, 'Position', [x 0.05 w h])
p3 = plot(dat.hrnmh, 100*dat.hrnmp, 'k.');
ylabel('detected (%)', 'Position', [8e-5, 50])
xscale log; xlabel('$h$','interpreter','latex', 'Position', [3e-3, -20]);
xlim([1e-4, 1e-1]); ylim([-10, 110]);  yticks([0, 100]);
text(0.0,1.1, '(c)','Units','normalized', 'FontSize', 12)

set(gcf,'Position',[3000 100 350 550])