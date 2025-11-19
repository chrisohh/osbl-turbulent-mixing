%% Power spectral density
%% pwelch 
% Take steady-state portion only
% time_raw=data_collection(6).time;
% u_mag=sqrt((data_collection(6).U).^2+(data_collection(6).V).^2+(data_collection(6).W).^2);
% u_raw=data_collection(6).T;%_mag;
fs = 1/mean(diff(time_raw),'omitnan');
figure;plot(time_raw,u_mag)
[~,start_idx]=min(abs(time_raw-0));
[~,end_idx]=min(abs(time_raw-60));
steady_state_idx=[start_idx:end_idx];
U_steady=u_raw(start_idx:end_idx);
time_steady = time_raw(steady_state_idx);

% Power spectral density
[psd, freq] = pwelch(U_steady - mean(U_steady), [], [], [], fs);

figure;
loglog(freq, psd, 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('PSD (m^2/s^2/Hz)');
title('Power Spectral Density - Steady State');


%% 
% % Turbulence typically shows:
% % - Low frequency content (large eddies): < 10 Hz for typical flows
% % - Inertial subrange: f^(-5/3) slope in the middle
% % - White noise plateau at high frequencies
% % 
% % Find where spectrum flattens (noise floor)
% log_freq = log10(freq_welch(freq_welch > 1));
% log_psd = log10(psd_welch(freq_welch > 1));
% 
% % Calculate local slope
% slope = diff(log_psd) ./ diff(log_freq);
% smooth_slope = movmean(slope, 10);
% 
% % Find where slope > -0.5 (approaching flat/noise region)
% noise_start_idx = find(smooth_slope > -0.5, 1);
% if ~isempty(noise_start_idx)
%     noise_freq = freq_welch(find(freq_welch > 1, 1) + noise_start_idx);
%     fprintf('Estimated noise frequency cutoff: %.1f Hz\n', noise_freq);
%     hold on;
%     xline(noise_freq, 'r--', 'Noise floor', 'LineWidth', 2);
% end

% % Apply low-pass filter to remove high-frequency noise
% cutoff_freq = 50; % Hz - adjust based on your spectrum
% [b, a] = butter(4, cutoff_freq/(fs/2), 'low');
% U_filtered = filtfilt(b, a, U_steady);

%% Simple moving average filter (acts as low-pass)
filter_window = round(fs / 460); % Cutoff ~460 Hz
U_filtered = movmean(U_steady, filter_window);

figure;
plot(time_steady, U_steady, 'Color', [0.7 0.7 0.7]);
hold on;
plot(time_steady, U_filtered, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Raw (noise + turbulence)', 'Filtered (turbulence only)');
title('Effect of Moving Average Filter');

% Calculate RMS of turbulence vs noise
turb_rms = std(U_filtered - mean(U_filtered));
noise_rms = std(U_steady - U_filtered);
total_rms = std(U_steady - mean(U_steady));

fprintf('\nFluctuation breakdown:\n');
fprintf('Total RMS: %.4f m/s\n', total_rms);
fprintf('Turbulent RMS: %.4f m/s (%.1f%%)\n', turb_rms, 100*turb_rms^2/total_rms^2);
fprintf('Noise RMS: %.4f m/s (%.1f%%)\n', noise_rms, 100*noise_rms^2/total_rms^2);

%% Reynolds Stress Analysis
% n = 1; % Select dataset

% % Use only steady-state data
% U_steady = data_collection(n).U(steady_state_idx);
% V_steady = data_collection(n).V(steady_state_idx);
% W_steady = data_collection(n).W(steady_state_idx);

% Optional: Apply low-pass filter to remove noise (recommended)
cutoff_freq = 10; % Hz - based on your PSD analysis
filter_window = round(fs / cutoff_freq);
U_steady = movmean(U_steady, filter_window);
V_steady = movmean(V_steady, filter_window);
W_steady = movmean(W_steady, filter_window);

% Calculate mean velocities
U_mean = mean(U_steady);
V_mean = mean(V_steady);
W_mean = mean(W_steady);

% Calculate fluctuating components (Reynolds decomposition)
u_prime = U_steady - U_mean;
v_prime = V_steady - V_mean;
w_prime = W_steady - W_mean;

% Normal stresses (variance of each component)
uu_stress = mean(u_prime.^2);  % <u'u'> or <u'^2>
vv_stress = mean(v_prime.^2);  % <v'v'> or <v'^2>
ww_stress = mean(w_prime.^2);  % <w'w'> or <w'^2>

% Shear stresses (covariances)
uv_stress = mean(u_prime .* v_prime);  % <u'v'>
uw_stress = mean(u_prime .* w_prime);  % <u'w'>
vw_stress = mean(v_prime .* w_prime);  % <v'w'>

% Turbulent kinetic energy
TKE = 0.5 * (uu_stress + vv_stress + ww_stress);

fprintf('\n=== Reynolds Stress Tensor ===\n');
fprintf('Normal stresses:\n');
fprintf('  <u''²> = %.6f m²/s²\n', uu_stress);
fprintf('  <v''²> = %.6f m²/s²\n', vv_stress);
fprintf('  <w''²> = %.6f m²/s²\n', ww_stress);
fprintf('\nShear stresses:\n');
fprintf('  <u''v''> = %.6f m²/s²\n', uv_stress);
fprintf('  <u''w''> = %.6f m²/s²\n', uw_stress);
fprintf('  <v''w''> = %.6f m²/s²\n', vw_stress);
fprintf('\nTurbulent Kinetic Energy (TKE) = %.6f m²/s²\n', TKE);

% RMS values (standard deviations)
u_rms = sqrt(uu_stress);
v_rms = sqrt(vv_stress);
w_rms = sqrt(ww_stress);

fprintf('\nRMS fluctuations:\n');
fprintf('  u_rms = %.4f m/s\n', u_rms);
fprintf('  v_rms = %.4f m/s\n', v_rms);
fprintf('  w_rms = %.4f m/s\n', w_rms);

% Turbulence intensities (percentage)
U_mag_mean = sqrt(U_mean^2 + V_mean^2 + W_mean^2);
TI_u = 100 * u_rms / U_mag_mean;
TI_v = 100 * v_rms / U_mag_mean;
TI_w = 100 * w_rms / U_mag_mean;
TI_total = 100 * sqrt((uu_stress + vv_stress + ww_stress)/3) / U_mag_mean;

fprintf('\nTurbulence intensities:\n');
fprintf('  TI_u = %.2f%%\n', TI_u);
fprintf('  TI_v = %.2f%%\n', TI_v);
fprintf('  TI_w = %.2f%%\n', TI_w);
fprintf('  TI_total = %.2f%%\n', TI_total);

%% Full Reynolds stress tensor
% τ_ij = ρ * <u_i' u_j'>
% For incompressible flow, often shown normalized by density or as specific

Reynolds_stress_tensor = [uu_stress, uv_stress, uw_stress;
                          uv_stress, vv_stress, vw_stress;
                          uw_stress, vw_stress, ww_stress];

fprintf('\nReynolds Stress Tensor:\n');
fprintf('         U          V          W\n');
fprintf('U  %10.6f %10.6f %10.6f\n', Reynolds_stress_tensor(1,:));
fprintf('V  %10.6f %10.6f %10.6f\n', Reynolds_stress_tensor(2,:));
fprintf('W  %10.6f %10.6f %10.6f\n', Reynolds_stress_tensor(3,:));

% Check if tensor is symmetric (it should be)
if max(max(abs(Reynolds_stress_tensor - Reynolds_stress_tensor'))) < 1e-10
    fprintf('\n✓ Tensor is symmetric (as expected)\n');
end

%% Anisotropy tensor (deviation from isotropic turbulence)
% b_ij = <u_i' u_j'> / (2*TKE) - δ_ij/3

anisotropy_tensor = Reynolds_stress_tensor / (2*TKE) - eye(3)/3;

fprintf('\nAnisotropy Tensor:\n');
fprintf('         U          V          W\n');
fprintf('U  %10.6f %10.6f %10.6f\n', anisotropy_tensor(1,:));
fprintf('V  %10.6f %10.6f %10.6f\n', anisotropy_tensor(2,:));
fprintf('W  %10.6f %10.6f %10.6f\n', anisotropy_tensor(3,:));

% Anisotropy invariants
II = -0.5 * trace(anisotropy_tensor^2);
III = det(anisotropy_tensor);

fprintf('\nAnisotropy invariants:\n');
fprintf('  II = %.6f\n', II);
fprintf('  III = %.6f\n', III);

% Perfect isotropic turbulence: II = 0, III = 0
% One-component turbulence: II = 2/9, III = 2/27

%% 1. Reynolds stress bar chart
figure;
stress_values = [uu_stress, vv_stress, ww_stress, uv_stress, uw_stress, vw_stress];
stress_labels = {'<u''²>', '<v''²>', '<w''²>', '<u''v''>', '<u''w''>', '<v''w''>'};
bar(stress_values);
set(gca, 'XTickLabel', stress_labels);
ylabel('Stress (m²/s²)');
title('Reynolds Stresses');
grid on;

%% 2. Time history of fluctuations
figure;
time_plot = data_collection(n).time(steady_state_idx);
subplot(3,1,1);
plot(time_plot, u_prime);
ylabel('u'' (m/s)'); title('Velocity Fluctuations');
grid on;

subplot(3,1,2);
plot(time_plot, v_prime);
ylabel('v'' (m/s)');
grid on;

subplot(3,1,3);
plot(time_plot, w_prime);
ylabel('w'' (m/s)');
xlabel('Time (s)');
grid on;

%% 3. PDF of fluctuations
figure;
subplot(1,3,1);
histogram(u_prime, 50, 'Normalization', 'pdf');
xlabel('u'' (m/s)'); ylabel('PDF');
title('u'' distribution');
grid on;

subplot(1,3,2);
histogram(v_prime, 50, 'Normalization', 'pdf');
xlabel('v'' (m/s)'); ylabel('PDF');
title('v'' distribution');
grid on;

subplot(1,3,3);
histogram(w_prime, 50, 'Normalization', 'pdf');
xlabel('w'' (m/s)'); ylabel('PDF');
title('w'' distribution');
grid on;

%% 4. Scatter plots of correlations
figure;
subplot(1,3,1);
scatter(u_prime, v_prime, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('u'' (m/s)'); ylabel('v'' (m/s)');
title(sprintf('<u''v''> = %.6f', uv_stress));
axis equal; grid on;

subplot(1,3,2);
scatter(u_prime, w_prime, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('u'' (m/s)'); ylabel('w'' (m/s)');
title(sprintf('<u''w''> = %.6f', uw_stress));
axis equal; grid on;

subplot(1,3,3);
scatter(v_prime, w_prime, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('v'' (m/s)'); ylabel('w'' (m/s)');
title(sprintf('<v''w''> = %.6f', vw_stress));
axis equal; grid on;

%% 5. Anisotropy invariant map (Lumley triangle)
figure;
hold on;

% Plot Lumley triangle boundaries
III_line = linspace(-1/27, 1/9, 100);
II_upper = (-3*III_line).^(1/2);
II_lower = (4/3 * III_line).^(1/3);

plot(III_line, II_upper, 'k-', 'LineWidth', 2);
plot(III_line, II_lower, 'k-', 'LineWidth', 2);

% Mark special points
plot(0, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % Isotropic
text(0, 0.02, 'Isotropic', 'HorizontalAlignment', 'center');

plot(2/27, 2/9, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % 1D
text(2/27, 0.24, '1D', 'HorizontalAlignment', 'center');

plot(-1/27, 1/9, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % 2D
text(-1/27, 0.14, '2D', 'HorizontalAlignment', 'center');

% Plot your data
plot(III, II, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'k');

xlabel('III'); ylabel('II');
title('Anisotropy Invariant Map (Lumley Triangle)');
grid on; axis equal;
xlim([-0.06, 0.1]); ylim([0, 0.25]);


%% Primary Reynolds stress for ocean boundary layer: <u'w'>
% This represents vertical momentum flux - key for air-sea interaction

figure('Position', [100, 100, 1200, 400]);

% Time series of instantaneous momentum flux
subplot(1,3,1);
time_plot = data_collection(n).time(steady_state_idx);
momentum_flux_instant = u_prime .* w_prime;
plot(time_plot, momentum_flux_instant, 'Color', [0.5 0.5 0.5]);
hold on;
plot(time_plot, movmean(momentum_flux_instant, round(fs*0.1)), 'b-', 'LineWidth', 2);
yline(uw_stress, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('u''w'' (m²/s²)');
title('Instantaneous Vertical Momentum Flux');
legend('Instantaneous', '0.1s average', 'Mean <u''w''>');
grid on;

% Quadrant analysis
subplot(1,3,2);
scatter(u_prime, w_prime, 10, momentum_flux_instant, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
xline(0, 'k-', 'LineWidth', 1);
yline(0, 'k-', 'LineWidth', 1);
xlabel('u'' (m/s)'); ylabel('w'' (m/s)');
title(sprintf('Momentum Flux: <u''w''> = %.6f m²/s²', uw_stress));
colorbar; colormap(jet);
axis equal;

% Add quadrant labels
xlims = xlim; ylims = ylim;
text(0.7*xlims(2), 0.7*ylims(2), 'Q1: Outward', 'FontSize', 10, 'Color', 'k');
text(0.7*xlims(1), 0.7*ylims(2), 'Q2: Ejection', 'FontSize', 10, 'Color', 'k');
text(0.7*xlims(1), 0.7*ylims(1), 'Q3: Inward', 'FontSize', 10, 'Color', 'k');
text(0.7*xlims(2), 0.7*ylims(1), 'Q4: Sweep', 'FontSize', 10, 'Color', 'k');
grid on;

% Quadrant contributions
subplot(1,3,3);
Q1 = sum(momentum_flux_instant(u_prime > 0 & w_prime > 0));
Q2 = sum(momentum_flux_instant(u_prime < 0 & w_prime > 0));
Q3 = sum(momentum_flux_instant(u_prime < 0 & w_prime < 0));
Q4 = sum(momentum_flux_instant(u_prime > 0 & w_prime < 0));
total = Q1 + Q2 + Q3 + Q4;

quadrant_contrib = [Q1, Q2, Q3, Q4] / total * 100;
bar(quadrant_contrib);
set(gca, 'XTickLabel', {'Q1: Outward', 'Q2: Ejection', 'Q3: Inward', 'Q4: Sweep'});
ylabel('Contribution to <u''w''> (%)');
title('Quadrant Analysis');
grid on;
xtickangle(45);

fprintf('\nQuadrant Analysis:\n');
fprintf('Q1 (Outward, u''>0, w''>0): %.1f%%\n', quadrant_contrib(1));
fprintf('Q2 (Ejection, u''<0, w''>0): %.1f%%\n', quadrant_contrib(2));
fprintf('Q3 (Inward, u''<0, w''<0): %.1f%%\n', quadrant_contrib(3));
fprintf('Q4 (Sweep, u''>0, w''<0): %.1f%%\n', quadrant_contrib(4));

%% Friction velocity (u_*) - critical for ocean surface layer
% u_* = sqrt(-<u'w'>) for downward momentum flux

if uw_stress < 0
    u_star = sqrt(-uw_stress);
    fprintf('\nFriction velocity u_* = %.4f m/s\n', u_star);
    
    % Surface stress (τ = ρ * u_*²)
    rho_air = 1.225; % kg/m³ at sea level, adjust for conditions
    tau_surface = rho_air * u_star^2;
    fprintf('Surface stress τ = %.6f Pa\n', tau_surface);
else
    fprintf('\nWarning: <u''w''> is positive - unexpected for downward momentum flux\n');
    u_star = sqrt(abs(uw_stress));
end

% Normalized stresses (by friction velocity)
uu_plus = uu_stress / u_star^2;
vv_plus = vv_stress / u_star^2;
ww_plus = ww_stress / u_star^2;

fprintf('\nNormalized Reynolds stresses:\n');
fprintf('<u''²>/u_*² = %.3f\n', uu_plus);
fprintf('<v''²>/u_*² = %.3f\n', vv_plus);
fprintf('<w''²>/u_*² = %.3f\n', ww_plus);

%% TKE production and dissipation estimates
% Production: P = -<u'w'> * dU/dz (for boundary layer)

figure;
subplot(2,2,1);
% TKE time series
TKE_instant = 0.5 * (u_prime.^2 + v_prime.^2 + w_prime.^2);
plot(time_plot, TKE_instant);
hold on;
yline(TKE, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('TKE (m²/s²)');
title('Turbulent Kinetic Energy');
legend('Instantaneous', 'Mean');
grid on;

subplot(2,2,2);
% Energy distribution
energy_components = [uu_stress, vv_stress, ww_stress] / (2*TKE) * 100;
bar(energy_components);
set(gca, 'XTickLabel', {'<u''²>', '<v''²>', '<w''²>'});
ylabel('% of Total TKE');
title('TKE Component Distribution');
grid on;

subplot(2,2,3);
% Anisotropy ratios (important for ocean surface layer)
ratio_vv_ww = vv_stress / ww_stress;
ratio_uu_ww = uu_stress / ww_stress;
bar([ratio_uu_ww, ratio_vv_ww, 1]);
set(gca, 'XTickLabel', {'<u''²>/<w''²>', '<v''²>/<w''²>', '<w''²>/<w''²>'});
ylabel('Stress Ratio');
title('Anisotropy Ratios');
yline(1, 'k--');
grid on;

fprintf('\nAnisotropy ratios:\n');
fprintf('<u''²>/<w''²> = %.3f\n', ratio_uu_ww);
fprintf('<v''²>/<w''²> = %.3f\n', ratio_vv_ww);

subplot(2,2,4);
% Structure parameter (A1 = 1 - 9/8*(uu-vv)²/TKE²)
A1 = 1 - 9/8 * (uu_stress - vv_stress)^2 / (2*TKE)^2;
bar([A1]);
ylim([0 1.2]);
set(gca, 'XTickLabel', {'A_1'});
ylabel('Value');
title(sprintf('Axisymmetric Structure Parameter\nA_1 = %.3f', A1));
text(1, A1+0.1, sprintf('A_1=%.3f\n(1=isotropic, 0=axisymmetric)', A1), ...
    'HorizontalAlignment', 'center');
grid on;