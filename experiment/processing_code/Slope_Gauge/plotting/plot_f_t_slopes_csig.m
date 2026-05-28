function plot_f_t_slopes_csig(data, frame_indices)
% PLOT_F_T_SLOPES_CSIG  Frequency-vs-time spectrogram of Sx and Sy.
%
% Same row-averaged spatial FFT as plot_kx_t_slopes_csig, then the
% one-sided kx axis is mapped to frequency via the linear gravity-capillary
% dispersion relation  omega = sqrt(g*kx + gamma*kx^3),  f = omega/(2*pi).
% The result is interpolated onto a uniform frequency grid for display.
%
% Also produces:
%   - peak-frequency ridge overlaid on the spectrogram
%   - separate figure: peak frequency f_peak(t) and spectral power at the
%     peak S(f_peak, t) for both Sx and Sy

Sx   = data.Sx;
Sy   = data.Sy;
time = data.time;
dx   = data.dx;

[Ny, Nx, ~] = size(Sx);

if nargin < 2 || isempty(frame_indices)
    idx = 1:size(Sx, 3);
else
    idx = frame_indices(:).';
end
n_t = numel(idx);

% One-sided kx axis (rad/m); kx(1)=0 is DC, skip when converting to f
Nk    = floor(Nx/2) + 1;
kx    = 2*pi * (0:Nk-1) / (Nx*dx);
win_x = hann(Nx).';

% Compute kx-t power spectra (Nk x n_t)
Sx_kt = zeros(Nk, n_t);
Sy_kt = zeros(Nk, n_t);
for j = 1:n_t
    [Sx_kt(:,j), Sy_kt(:,j)] = row_avg_spectrum( ...
        double(Sx(:,:,idx(j))), double(Sy(:,:,idx(j))), win_x, Nk);
end

t_plot = time(idx);

% --- Convert kx axis to frequency (Hz) via gravity-capillary dispersion ---
g     = 9.81;
gamma = 7.4e-5;          % sigma/rho for water, m^3 s^-2
k_use = kx(2:end);       % drop DC
f_kx  = sqrt(g .* k_use + gamma .* k_use.^3) / (2*pi);  % Hz, non-uniform

% Uniform f grid spanning same range
Nf    = Nk - 1;
f_uni = linspace(f_kx(1), f_kx(end), Nf);

% Interpolate each column (time snapshot) onto uniform f grid
Sx_ft = interp1(f_kx(:), Sx_kt(2:end,:), f_uni(:), 'linear', 0);  % Nf x n_t
Sy_ft = interp1(f_kx(:), Sy_kt(2:end,:), f_uni(:), 'linear', 0);

% --- Peak frequency: search above 0.5 Hz to avoid DC leakage ---
f_lo  = 0.5;
i_lo  = find(f_uni >= f_lo, 1);

[peak_pow_Sx, peak_idx_Sx] = max(Sx_ft(i_lo:end, :), [], 1);
[peak_pow_Sy, peak_idx_Sy] = max(Sy_ft(i_lo:end, :), [], 1);
peak_f_Sx = f_uni(peak_idx_Sx + i_lo - 1);
peak_f_Sy = f_uni(peak_idx_Sy + i_lo - 1);

%% f-t spectrogram
fig = figure('Name', 'f-t spectrogram (S_x, S_y)', ...
             'Position', [100 100 1100 800], 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
imagesc(t_plot, f_uni, log10(Sx_ft + eps));
set(gca, 'YDir', 'normal'); axis tight;
caxis([0, 4]);
colormap(ax1, bone);
cb = colorbar;
cb.Label.String      = '$\log_{10}|S_x(f)|^2$';
cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;
hold on;
plot(t_plot, peak_f_Sx, 'r-', 'LineWidth', 1.5);
xlabel('$t$ (s)',    'Interpreter', 'latex');
ylabel('$f$ (Hz)',   'Interpreter', 'latex');
title('$|S_x(f,t)|^2$ ($y$-averaged, gravity-capillary dispersion)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');

ax2 = nexttile;
imagesc(t_plot, f_uni, log10(Sy_ft + eps));
set(gca, 'YDir', 'normal'); axis tight;
caxis([0, 4]);
colormap(ax2, bone);
cb = colorbar;
cb.Label.String      = '$\log_{10}|S_y(f)|^2$';
cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;
hold on;
plot(t_plot, peak_f_Sy, 'r-', 'LineWidth', 1.5);
xlabel('$t$ (s)',    'Interpreter', 'latex');
ylabel('$f$ (Hz)',   'Interpreter', 'latex');
title('$|S_y(f,t)|^2$ ($y$-averaged, gravity-capillary dispersion)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');

%% Peak frequency and peak power evolution
fig2 = figure('Name', 'Peak frequency and power evolution', ...
              'Position', [100 100 1000 600], 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t_plot, peak_f_Sx, 'b-', 'LineWidth', 1.5); hold on;
plot(t_plot, peak_f_Sy, 'r-', 'LineWidth', 1.5);
xlabel('$t$ (s)',                  'Interpreter', 'latex');
ylabel('$f_\mathrm{peak}$ (Hz)',   'Interpreter', 'latex');
title('Peak frequency vs time',    'Interpreter', 'latex');
legend('$S_x$', '$S_y$',          'Interpreter', 'latex', 'Location', 'best');
set(gca, 'fontsize', 14, 'fontname', 'times');

nexttile;
plot(t_plot, peak_pow_Sx, 'b-', 'LineWidth', 1.5); hold on;
plot(t_plot, peak_pow_Sy, 'r-', 'LineWidth', 1.5);
xlabel('$t$ (s)',                               'Interpreter', 'latex');
ylabel('$|S(f_\mathrm{peak})|^2$',             'Interpreter', 'latex');
title('Spectral peak power vs time',            'Interpreter', 'latex');
legend('$S_x$', '$S_y$',                       'Interpreter', 'latex', 'Location', 'best');
set(gca, 'fontsize', 14, 'fontname', 'times');

%% Sanity-check: per-row spectra in frequency domain at first snapshot
j_check = 1;
[~, ~, Px_rows, Py_rows] = row_avg_spectrum( ...
    double(Sx(:,:,idx(j_check))), double(Sy(:,:,idx(j_check))), win_x, Nk);

Px_rows_f = interp1(f_kx(:), Px_rows(:, 2:end).', f_uni(:), 'linear', 0);
Py_rows_f = interp1(f_kx(:), Py_rows(:, 2:end).', f_uni(:), 'linear', 0);

n_show = min(8, Ny);
y_show = round(linspace(1, Ny, n_show));

fig3 = figure('Name', 'f-domain per-row sanity check', ...
              'Position', [100 100 1200 500], 'Color', 'w');
subplot(1,2,1);
loglog(f_uni(2:end), Px_rows_f(2:end, y_show), 'Color', [.6 .6 .6]); hold on;
loglog(f_uni(2:end), mean(Px_rows_f(2:end,:), 2), 'k', 'LineWidth', 2);
xlabel('f (Hz)'); ylabel('|S_x(f)|^2');
title(sprintf('S_x per-row (grey) vs y-average (black) at t = %.2f s', t_plot(j_check)));
grid on;

subplot(1,2,2);
loglog(f_uni(2:end), Py_rows_f(2:end, y_show), 'Color', [.6 .6 .6]); hold on;
loglog(f_uni(2:end), mean(Py_rows_f(2:end,:), 2), 'k', 'LineWidth', 2);
xlabel('f (Hz)'); ylabel('|S_y(f)|^2');
title(sprintf('S_y per-row (grey) vs y-average (black) at t = %.2f s', t_plot(j_check)));
grid on;
end


function [px_mean, py_mean, Px_rows, Py_rows] = row_avg_spectrum(fx, fy, win_x, Nk)
fx_d    = fx - mean(fx, 2);
fy_d    = fy - mean(fy, 2);
Fx      = fft(fx_d .* win_x, [], 2);
Fy      = fft(fy_d .* win_x, [], 2);
Px_rows = abs(Fx(:, 1:Nk)).^2;
Py_rows = abs(Fy(:, 1:Nk)).^2;
px_mean = mean(Px_rows, 1).';
py_mean = mean(Py_rows, 1).';
end
