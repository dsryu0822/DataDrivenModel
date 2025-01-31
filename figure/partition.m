dat = readtable("G:/DDM/partition.csv");

subplot(1,3,1)
p1 = plot(dat.softh, 100*dat.softp, 'k.');
ylabel('? (%)')
xscale log; xlabel('$h$','interpreter','latex');
ylim([-1, 101]); yticks([0, 50, 100]);
text(0.0,1.05, '(a)','Units','normalized')

subplot(1,3,2)
p2 = plot(dat.gearh, 100*dat.gearp, 'k.');
xscale log; xlabel('$h$','interpreter','latex');
ylim([-1, 101]); yticks([0, 50, 100]);
text(0.0,1.05, '(b)','Units','normalized')

subplot(1,3,3)
p3 = plot(dat.hrnmh, 100*dat.hrnmp, 'k.');
xscale log; xlabel('$h$','interpreter','latex');
ylim([-1, 101]); yticks([0, 50, 100]);
text(0.0,1.05, '(c)','Units','normalized')

set(gcf,'Position',[100 100 1200 300])