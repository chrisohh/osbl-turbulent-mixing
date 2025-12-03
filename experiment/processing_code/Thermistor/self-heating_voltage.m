%% Bridge Voltage vs Self-Heating for Various Thermistor Resistances
clear; clc; close all;

%% Parameters
R0_values = [2, 8, 10, 20] * 1000;  % Thermistor resistances (Ω)
DC_water = 0.25e-3;                        % Dissipation constant (W/°C)

%% Bridge voltage range
E_bridge = linspace(0.05, 1.5, 200);  % 0.05 to 1.5V

%% Calculate self-heating for each configuration
n_configs = length(R0_values);
self_heating = zeros(n_configs, length(E_bridge));

for i = 1:n_configs
    RT = R0_values(i);           % Thermistor resistance at 25°C
    R_bridge = RT * 1.1;         % Bridge resistance (roughly matched)
    
    for j = 1:length(E_bridge)
        E = E_bridge(j);
        
        % Current through thermistor
        I = E / (RT + R_bridge);
        
        % Power dissipated
        P = I^2 * RT;
        
        % Self-heating (in mK)
        self_heating(i, j) = (P / DC_water) * 1e3;  % mK
    end
end

%% Create plot
figure('Position', [100, 100, 1000, 700]);

% Plot curves
colors = lines(n_configs);
hold on;

for i = 1:n_configs
    plot(E_bridge, self_heating(i,:), 'LineWidth', 2.5, 'Color', colors(i,:));
end

% Reference lines
plot([0.01,10],[5,5],'k--','linewidth',1.5)
plot([0.01,10],[50,50],'k--','linewidth',1.5)
plot([0.01,10],[500,500],'k--','linewidth',1.5)


% Labels
set(gca,'fontsize',14)
xlabel('Bridge Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Self-heating (mK)', 'FontSize', 14, 'FontWeight', 'bold');
title('Self-Heating vs Bridge Voltage @ 25°C', 'FontSize', 15, 'FontWeight', 'bold');

% Legend
legend_labels = arrayfun(@(r) sprintf('%gkΩ @ 25°C', r/1000), R0_values, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northwest', 'FontSize', 11);

% Annotation

text(0.6, 6, '< 5 mK for 0.01^\circC accuracy', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
text(0.6, 60, '< 50 mK for 0.1^\circC accuracy', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
text(0.6, 600, '< 500 mK for 1^\circC accuracy', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');

grid on;
set(gca, 'xscale','log','yscale','log')
xlim([0.1,1])
ylim([0.5,500])
set(gca,'xtick',[0.1:0.1:1])
hold off;
%% tolerance
T = 25;                                    % Temperature (°C)
for R0 = R0_values
    RT = thermistor_resistance(25, R0, beta, T0);
    
    % dR/dT (negative for NTC)
    dR_dT = -RT * beta / (T+273.15)^2;
    
    % Fractional sensitivity
    fractional = abs(dR_dT) / RT;
    
    fprintf('%6.0f kΩ   | %8.0f | %8.1f | %9.4f | %.4f %%/°C\n', ...
            R0/1000, RT, dR_dT, fractional, fractional*100);
end