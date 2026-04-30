%% FP07 Thermistor Resistance vs Temperature Plot
clear; clc; close all;

%% Parameters
R0 = [2,8,10,20]*1000;     % 10kΩ at 25°C
beta = 3521;    % Material constant (K) B25/85
T0 = 25;        % Reference temperature (°C)

%% Temperature range for plotting
T_range = linspace(0, 40, 1000);  % 0 to 40°C, 1000 points

%% Calculate resistance curve
R_range = thermistor_resistance(T_range, R0', beta, T0);

%% Your specific operating points
T_operating = [18, 30];
R_operating = thermistor_resistance(T_operating, R0', beta, T0);

%% Plot
figure('Position', [100, 100, 800, 600]);

% Main curve
plot(T_range, R_range/1000, 'LineWidth', 2);
hold on;

% Operating points
plot(T_operating, R_operating/1000, 'ko', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'k', 'LineWidth', 2);

% % Highlight your water temperature range
% patch([18 20 20 18], [0 0 15 15], 'g', ...
%       'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Labels and formatting
set(gca,'fontsize',14)
xlabel('Temperature (°C)', 'FontSize', 12);
ylabel('Resistance (k$\Omega$)', 'FontSize', 14,'Interpreter','latex');
title('Various Nominal Resistances of FP07 Thermistors', 'FontSize', 14);
grid on;
legend_labels = arrayfun(@(r) sprintf('%gkΩ @ 25°C', r/1000), R0_values, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northwest', 'FontSize', 11);

% Add text annotations for operating points
for i = 1:size(R_operating,2)
    for j=1:size(R_operating,1)
    text(T_operating(i)+2.5, R_operating(j,i)/1000, ...
         sprintf('%.1f k\\Omega', R_operating(j,i)/1000), ...
         'FontSize', 12, 'HorizontalAlignment', 'center');
    end
end

hold off;