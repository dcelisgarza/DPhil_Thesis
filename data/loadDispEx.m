close all
%% 100
Ni100 = readtable('Ni100.csv', 'PreserveVariableNames', false);
Ni100 = table2array(Ni100(2:end, 1:5));
Ni100_5a = rmmissing(Ni100(:, 1:2));
Ni100_5b = rmmissing(Ni100(:, [1, 3]));
Ni100_500a = rmmissing(Ni100(:, [1, 4]));
Ni100_500b = rmmissing(Ni100(:, [1, 5]));

%% 110
Ni110 = readtable('Ni110.csv', 'PreserveVariableNames', false);
Ni110 = table2array(Ni110(2:end, 1:5));
Ni110_5a = rmmissing(Ni110(:, 1:2));
Ni110_5b = rmmissing(Ni110(:, [1, 3]));
Ni110_500a = rmmissing(Ni110(:, [1, 4]));
Ni110_500b = rmmissing(Ni110(:, [1, 5]));

%% 100
fig1 = figure;
hold on
plot(Ni100_5a(:, 1), Ni100_5a(:, 2), '-', 'LineWidth', 2)
plot(Ni100_5b(:, 1), Ni100_5b(:, 2), '--', 'LineWidth', 2)
plot(Ni100_500a(:, 1), Ni100_500a(:, 2), '-.', 'LineWidth', 2)
plot(Ni100_500b(:, 1), Ni100_500b(:, 2), ':', 'LineWidth', 2)
xlabel('Strain, \mum/\mum', 'FontSize', 18)
ylabel('Stress, MPa', 'FontSize', 18)
legend('(a) 5 nm/s', '(b) 5 nm/s', '(a) 500 nm/s', '(b) 500 nm/s', 'FontSize', 15, 'Location', 'best')
ax = gca;
ax.FontSize = 15;
set(gca,'YLim',[0 300],'YTick',[0:25:300])
hold off
axis tight manual
set(fig1, 'Units', 'Inches');
pos = get(fig1, 'Position');
set(fig1, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig1, 'Ni100.pdf', '-dpdf', '-r300')

%% 110
fig2 = figure;
hold on
plot(Ni110_5a(:, 1), Ni110_5a(:, 2), '-', 'LineWidth', 2)
plot(Ni110_5b(:, 1), Ni110_5b(:, 2), '--', 'LineWidth', 2)
plot(Ni110_500a(:, 1), Ni110_500a(:, 2), '-.', 'LineWidth', 2)
plot(Ni110_500b(:, 1), Ni110_500b(:, 2), ':', 'LineWidth', 2)
xlabel('Strain, \mum/\mum', 'FontSize', 18)
ylabel('Stress, MPa', 'FontSize', 18)
legend('(a) 5 nm/s', '(b) 5 nm/s', '(a) 500 nm/s', '(b) 500 nm/s', 'FontSize', 15, 'Location', 'best')
ax = gca;
ax.FontSize = 15;
set(gca,'YLim',[0 200],'YTick',[0:20:200])
hold off
axis tight manual
set(fig2, 'Units', 'Inches');
pos = get(fig2, 'Position');
set(fig2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(fig2, 'Ni110.pdf', '-dpdf', '-r300')
