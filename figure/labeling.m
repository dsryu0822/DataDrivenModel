dat = readtable("G:/DDM/labeling.csv");

plot(dat.x_, 100*dat.acc, 'k.')
xscale log; xlabel('$\theta$','interpreter','latex');
xticks([1e-20, 1e+2]); xlim([1e-20, 1e+2]);
ylabel('accuracy (%)')
yticks([99, 100]); ylim([98.8, 100.1]);