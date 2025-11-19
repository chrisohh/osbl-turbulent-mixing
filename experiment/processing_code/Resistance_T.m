%% FP07 Thermistor Resistance vs Temperature Plot
clear; clc; close all;

%% Parameters
R0 = 10000;     % 10kΩ at 25°C
beta = 3550;    % Material constant (K)
T0 = 25;        % Reference temperature (°C)

%% Temperature range for plotting
T_range = linspace(0, 40, 1000);  % 0 to 40°C, 1000 points

%% Calculate resistance curve
R_range = thermistor_resistance(T_range, R0, beta, T0);

%% Your specific operating points
T_operating = [18, 19, 20, 25];
R_operating = thermistor_resistance(T_operating, R0, beta, T0);

%% Plot
figure('Position', [100, 100, 800, 600]);

% Main curve
plot(T_range, R_range/1000, 'b-', 'LineWidth', 2);
hold on;

% Operating points
plot(T_operating, R_operating/1000, 'ro', 'MarkerSize', 10, ...
     'MarkerFaceColor', 'r', 'LineWidth', 2);

% Highlight your water temperature range
patch([18 20 20 18], [0 0 15 15], 'g', ...
      'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Labels and formatting
xlabel('Temperature (°C)', 'FontSize', 12);
ylabel('Resistance (kΩ)', 'FontSize', 12);
title(sprintf('FP07 Thermistor: R_0 = %.0f Ω at %.0f°C, \\beta = %.0f K', ...
              R0, T0, beta), 'FontSize', 14);
grid on;
legend('R(T) Curve', 'Operating Points', 'Your Water Range (18-20°C)', ...
       'Location', 'northeast');

% Add text annotations for operating points
for i = 1:length(T_operating)
    text(T_operating(i), R_operating(i)/1000 + 0.3, ...
         sprintf('%.1f°C: %.2f kΩ', T_operating(i), R_operating(i)/1000), ...
         'FontSize', 9, 'HorizontalAlignment', 'center');
end

hold off;