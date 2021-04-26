close all
clear all
%%

% 100
Ni100 = readtable('Ni100.csv', 'PreserveVariableNames', false);
Ni100 = table2array(Ni100(2:end, 1:5));
Ni100_5a = rmmissing(Ni100(:, 1:2));
Ni100_5b = rmmissing(Ni100(:, [1, 3]));
Ni100_500a = rmmissing(Ni100(:, [1, 4]));
Ni100_500b = rmmissing(Ni100(:, [1, 5]));
[Ni100_pca,Ni100_score,Ni100_latent,Ni100_tsquared,Ni100_explained,Ni100_mu] = pca(Ni100(:, 2:end));

% 110
Ni110 = readtable('Ni110.csv', 'PreserveVariableNames', false);
Ni110 = table2array(Ni110(2:end, 1:5));
Ni110_5a = rmmissing(Ni110(:, 1:2));
Ni110_5b = rmmissing(Ni110(:, [1, 3]));
Ni110_500a = rmmissing(Ni110(:, [1, 4]));
Ni110_500b = rmmissing(Ni110(:, [1, 5]));
[Ni110_pca,Ni110_score,Ni110_latent,Ni110_tsquared,Ni110_explained,Ni110_mu] = pca(Ni110(:, 2:end));

%% 100
fig1 = figure;
hold on
plot(Ni100_5a(:, 1), Ni100_5a(:, 2), '-', 'LineWidth', 2)
plot(Ni100_5b(:, 1), Ni100_5b(:, 2), '--', 'LineWidth', 2)
plot(Ni100_500a(:, 1), Ni100_500a(:, 2), '-.', 'LineWidth', 2)
plot(Ni100_500b(:, 1), Ni100_500b(:, 2), ':', 'LineWidth', 2)
hold off
xlabel('Strain, \mum/\mum', 'FontSize', 18)
ylabel('Stress, MPa', 'FontSize', 18)
legend('(a) 5 nm/s', '(b) 5 nm/s', '(c) 500 nm/s', '(d) 500 nm/s', 'FontSize', 15, 'Location', 'best')
ax = gca;
ax.FontSize = 15;
set(gca,'YLim',[0 300],'YTick',[0:25:300])
axis tight manual
set(fig1, 'Units', 'Inches');
pos = get(fig1, 'Position');
set(fig1, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig1, 'Ni100.pdf', '-dpdf', '-r300')

% 110
fig2 = figure;
hold on
plot(Ni110_5a(:, 1), Ni110_5a(:, 2), '-', 'LineWidth', 2)
plot(Ni110_5b(:, 1), Ni110_5b(:, 2), '--', 'LineWidth', 2)
plot(Ni110_500a(:, 1), Ni110_500a(:, 2), '-.', 'LineWidth', 2)
plot(Ni110_500b(:, 1), Ni110_500b(:, 2), ':', 'LineWidth', 2)
xlabel('Strain, \mum/\mum', 'FontSize', 18)
ylabel('Stress, MPa', 'FontSize', 18)
legend('(a) 5 nm/s', '(b) 5 nm/s', '(c) 500 nm/s', '(d) 500 nm/s', 'FontSize', 15, 'Location', 'best')
ax = gca;
ax.FontSize = 15;
set(gca,'YLim',[0 200],'YTick',[0:20:200])
hold off
axis tight manual
set(fig2, 'Units', 'Inches');
pos = get(fig2, 'Position');
set(fig2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig2, 'Ni110.pdf', '-dpdf', '-r300')

%%
fig3 = figure;
hold on
scatter(Ni100_pca(1:2,1), Ni100_pca(1:2,2), 80, '^', 'filled')
scatter(Ni100_pca(3:4,1), Ni100_pca(3:4,2), 80, 'filled')
legend('(a) 5 nm/s', '(b) 500 nm/s', 'FontSize', 15, 'Location', 'best')
xlabel('PC 1', 'FontSize', 18)
ylabel('PC 2', 'FontSize', 18)
ax = gca;
ax.FontSize = 15;
hold off
% axis tight manual
set(fig3, 'Units', 'Inches');
pos = get(fig3, 'Position');
set(fig3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig3, 'Ni100_pca.pdf', '-dpdf', '-r300')

%%
fig4 = figure;
hold on
scatter(Ni110_pca(1:2,1), Ni110_pca(1:2,2), 80, '^', 'filled')
scatter(Ni110_pca(3:4,1), Ni110_pca(3:4,2), 80, 'filled')
legend('(a) 5 nm/s', '(b) 500 nm/s', 'FontSize', 15, 'Location', 'best')
xlabel('PC 1', 'FontSize', 18)
ylabel('PC 2', 'FontSize', 18)
ax = gca;
ax.FontSize = 15;
hold off
% axis tight manual
set(fig4, 'Units', 'Inches');
pos = get(fig4, 'Position');
set(fig4, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig4, 'Ni110_pca.pdf', '-dpdf', '-r300')

%%



addpath 'D:\DPhil\OneDrive - Nexus365\EasyDD\src'
Ni100a_ini = load('initial_11-Mar-2021_8_tensile_ni_100_0.mat');
Ni100a = load('11-Mar-2021_8_tensile_ni_100_143800.mat');
Ni100aN = load('11-Mar-2021_numT_8_tensile_ni_100_131400.mat');
Ni100b_ini = load('initial_16-Mar-2021_8_tensile_ni_100_0.mat');
Ni100b = load('16-Mar-2021_8_tensile_ni_100_141200.mat');

amag = Ni100a.amag;
mumag = Ni100a.mumag;
dx = Ni100a.dx;
dy = Ni100a.dy;
dz = Ni100a.dz;

Usim100a = Ni100a.Usim;
Fsim100a = Ni100a.Fsim;
curstep100a = Ni100a.curstep;

Usim100aN = Ni100aN.Usim;
Fsim100aN = Ni100aN.Fsim;
curstep100aN = Ni100aN.curstep;

Usim100b = Ni100b.Usim;
Fsim100b = Ni100b.Fsim;
curstep100b = Ni100b.curstep;

fig5 = figure;
hold on
plot(Usim100a(1:curstep100a)*amag/(dx*amag),Fsim100a(1:curstep100a)*(amag*1e-6)^2*mumag/(dz*amag*1e-6*dy*amag*1e-6), '-', 'LineWidth', 2)
plot(Usim100aN(1:curstep100aN)*amag/(dx*amag),Fsim100aN(1:curstep100aN)*(amag*1e-6)^2*mumag/(dz*amag*1e-6*dy*amag*1e-6), ':', 'LineWidth', 2)
plot(Usim100b(1:curstep100b)*amag/(dx*amag),Fsim100b(1:curstep100b)*(amag*1e-6)^2*mumag/(dz*amag*1e-6*dy*amag*1e-6), '-.', 'LineWidth', 2)
hold off
xlabel('Strain, \mum/\mum', 'FontSize', 18)
ylabel('Stress, MPa', 'FontSize', 18)
legend('(a) 2.5 cm/s', '(b) 2.5 cm/s', '(c) 2.5 cm/s', 'FontSize', 15, 'Location', 'best')
ax = gca;
ax.FontSize = 15;
set(gca,'YLim',[0 260],'YTick',[0:20:260])
set(gca,'XLim',[0 0.002],'XTick',[0:0.00025:0.002])
axis tight manual
set(fig5, 'Units', 'Inches');
pos = get(fig5, 'Position');
set(fig5, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig5, 'Ni100_DDD.pdf', '-dpdf', '-r300')


% [Ni100_pca,Ni100_score,Ni100_latent,Ni100_tsquared,Ni100_explained,Ni100_mu] = pca(Ni100(:, 2:end));
% 
% hi = min([max(Usim100a); max(Usim100aN)]*amag/(dx*amag));
% xq = linspace(0, hi, curstep100a);
% 
% F100a_i = interp1(Usim100a(1:curstep100a)*amag/(dx*amag),Fsim100a(1:curstep100a)*(amag*1e-6)^2*mumag/(dz*amag*1e-6*dy*amag*1e-6),xq);
% F100aN_i = interp1(Usim100aN(1:curstep100aN)*amag/(dx*amag),Fsim100aN(1:curstep100aN)*(amag*1e-6)^2*mumag/(dz*amag*1e-6*dy*amag*1e-6),xq);
% hold on
% plot(xq, F100a_i)
% plot(xq, F100aN_i)
% hold off
% 
% pca([F100a_i' F100aN_i'])
%%

Ni100b_ini = load('initial_16-Mar-2021_8_tensile_ni_100_0.mat');
Ni100b = load('16-Mar-2021_8_tensile_ni_100_141200.mat');

Ni110a_ini = load('initial_11-Mar-2021_4_tensile_ni_110_0.mat');
Ni110a = load('11-Mar-2021_4_tensile_ni_110_163800.mat');

Ni110b_ini = load('initial_16-Mar-2021_4_tensile_ni_110_0.mat');
Ni110b = load('16-Mar-2021_4_tensile_ni_110_154800.mat');

amag = Ni110a.amag;
mumag = Ni110a.mumag;
dx = Ni110a.dx;
dy = Ni110a.dy;
dz = Ni110a.dz;

Usim110a = Ni110a.Usim;
Fsim110a = Ni110a.Fsim;
curstep110a = Ni110a.curstep;

Usim110b = Ni110b.Usim;
Fsim110b = Ni110b.Fsim;
curstep110b = Ni110b.curstep;

fig6 = figure;
hold on
plot(Usim110a(1:curstep110a)*amag/(dx*amag),Fsim110a(1:curstep110a)*(amag*1e-6)^2*mumag/(dz*amag*1e-6*dy*amag*1e-6), '-', 'LineWidth', 2)
plot(Usim110b(1:curstep110b)*amag/(dx*amag),Fsim110b(1:curstep110b)*(amag*1e-6)^2*mumag/(dz*amag*1e-6*dy*amag*1e-6), '-.', 'LineWidth', 2)
hold off
xlabel('Strain, \mum/\mum', 'FontSize', 18)
ylabel('Stress, MPa', 'FontSize', 18)
legend('(a) 2.5 cm/s', '(b) 2.5 cm/s', 'FontSize', 15, 'Location', 'best')
ax = gca;
ax.FontSize = 15;
set(gca,'YLim',[0 260],'YTick',[0:20:260])
set(gca,'XLim',[0 0.002],'XTick',[0:0.00025:0.002])
axis tight manual
set(fig6, 'Units', 'Inches');
pos = get(fig6, 'Position');
set(fig6, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig6, 'Ni110_DDD.pdf', '-dpdf', '-r300')
%%
