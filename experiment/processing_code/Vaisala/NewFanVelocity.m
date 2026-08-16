Voltage = [2, 2.5, 3, 4, 5, 6, 7, 8, 8.5];
Velocity = [1.9, 2.3, 2.9, 4, 5.4, 6.5, 7.5, 8.9, 9.3];

% Linear regression
p = polyfit(Voltage, Velocity, 1);
Vfit = linspace(min(Voltage), max(Voltage), 100);
Velfit = polyval(p, Vfit);

% Plot
figure;
plot(Voltage, Velocity, 'o', 'markersize', 8, 'LineWidth', 1);
hold on;
plot(Vfit, Velfit, 'LineWidth', 1.5);

% Equation string
eqStr = sprintf('Velocity = %.2fV + %.2f', p(1), p(2));
text(2.2, 8.5, eqStr, 'FontSize', 13, 'FontName', 'times');
disp(eqStr)
xlabel('Voltage (V)');
ylabel('Velocity (m/s)');
% legend('Data', 'Linear fit', 'Location', 'northwest');
set(gca, 'FontSize', 14, 'FontName', 'times');
hold off;