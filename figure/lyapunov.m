%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               soft impact                     %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = readtable("G:/DDM/result/soft_lyapunov.csv");
subplot(2, 3, 1)
plot(dat.bp, dat.x_1, 'k-')
hold on
plot(dat.bp, dat.x_2, 'k-')
plot(dat.bp, dat.x_3, 'k-')
xticks([]);
yticks(0)
ylabel0 = ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel0.Position(1) = ylabel0.Position(1) + 0.012;
set(gca, 'Position', [0.03 0.51 0.3 0.45])
text(0.0,1.05, '(a)','Units','normalized', 'FontSize', 14)

dat = readtable("G:/DDM/result/soft_lyapunov_rcvd.csv");
subplot(2, 3, 4)
plot(dat.bp, dat.x_1, 'b-')
hold on
plot(dat.bp, dat.x_2, 'b-')
plot(dat.bp, dat.x_3, 'b-')
xticks([min(dat.bp), max(dat.bp)]);
yticks(0)
ylabel0 = ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel0.Position(1) = ylabel0.Position(1) + 0.012;
xlabel0 = xlabel('$d$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0.Position(2) = xlabel0.Position(2) + 0.3;
set(gca, 'Position', [0.03 0.05 0.3 0.45])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               gear system                     %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dat = readtable("G:/DDM/result/gear_lyapunov.csv");
subplot(2, 3, 2)
plot(dat.bp, dat.x_1, 'k-')
hold on
plot(dat.bp, dat.x_2, 'k-')
plot(dat.bp, dat.x_3, 'k-')
xticks([]);
yticks(0); % ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'Position', [0.36 0.51 0.3 0.45])
text(0.0,1.05, '(b)','Units','normalized', 'FontSize', 14)
ylim([-0.07, 0.03])

dat = readtable("G:/DDM/result/gear_lyapunov_rcvd.csv");
subplot(2, 3, 5)
plot(dat.bp, dat.x_1, 'b-')
hold on
plot(dat.bp, dat.x_2, 'b-')
plot(dat.bp, dat.x_3, 'b-')
xticks([min(dat.bp), max(dat.bp)]);
yticks(0); % ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0 = xlabel('$\bar{F}_{e}$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0.Position(2) = xlabel0.Position(2) - 0.002;
set(gca, 'Position', [0.36 0.05 0.3 0.45])
ylim([-0.07, 0.03])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%               Hindmarsh-Rose                  %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat = readtable("G:/DDM/result/hrnm_lyapunov.csv");
subplot(2, 3, 3)
plot(dat.bp, dat.x_1, 'k-')
hold on
plot(dat.bp, dat.x_2, 'k-')
plot(dat.bp, dat.x_3, 'k-')
xticks([]);
yticks(0); % ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'Position', [0.69 0.51 0.3 0.45])
text(0.0,1.05, '(c)','Units','normalized', 'FontSize', 14)

dat = readtable("G:/DDM/result/hrnm_lyapunov_rcvd.csv");
subplot(2, 3, 6)
plot(dat.bp, dat.x_1, 'b-')
hold on
plot(dat.bp, dat.x_2, 'b-')
plot(dat.bp, dat.x_3, 'b-')
xticks([min(dat.bp), max(dat.bp)]);
yticks(0); % ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0 = xlabel('$l$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel0.Position(2) = xlabel0.Position(2) - 0.01;
set(gca, 'Position', [0.69 0.05 0.3 0.45])

set(gcf,'Position', [3000 100 1000 600])
