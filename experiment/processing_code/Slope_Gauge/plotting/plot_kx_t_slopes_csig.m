function plot_kx_t_slopes_csig(data, frame_indices)
% PLOT_KX_T_SLOPES_CSIG  kx-vs-time spectrogram of Sx and Sy.
%
% For each snapshot t:
%   - For each y row: 1D FFT of Sx in the along-wind (x) direction with
%     a Hann window; same for Sy.
%   - Average the per-row |FFT|^2 across y -> one wavenumber spectrum.
% Stack across snapshots -> (kx, t) intensity map.
%
% Optional frame_indices: time indices into data.Sx (default = all).
% Also produces a sanity-check figure showing the per-row spectra and
% the y-average for the first plotted snapshot.

Sx   = data.Sx;
Sy   = data.Sy;
time = data.time;
dx   = data.dx;                          % m

[Ny, Nx, Nt] = size(Sx);

if nargin < 2 || isempty(frame_indices)
    idx = 1:Nt;
else
    idx = frame_indices(:).';
end
n_t = numel(idx);

% One-sided wavenumber axis (rad/m). DC at index 1, Nyquist at index Nk.
Nk = floor(Nx/2) + 1;
kx = 2*pi * (0:Nk-1) / (Nx*dx);

win_x = hann(Nx).';                      % row window applied along x

Sx_kt = zeros(Nk, n_t);
Sy_kt = zeros(Nk, n_t);

for j = 1:n_t
    [Sx_kt(:, j), Sy_kt(:, j)] = row_avg_spectrum( ...
        double(Sx(:,:,idx(j))), double(Sy(:,:,idx(j))), win_x, Nk);
end

t_plot = time(idx);

%% Main figure: kx–t intensity for Sx and Sy
fig = figure('Name', 'kx-t spectrogram (S_x, S_y)', ...
             'Position', [100 100 1100 800], 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
imagesc(t_plot, kx, log10(Sx_kt + eps));
set(gca, 'YDir', 'normal'); axis tight;
 caxis([0,4])
colormap(ax1, bone);
cb = colorbar; cb.Label.String = '$\log_{10}|S_x(k_x)|^2$';cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;

xlabel('$t$ (s)','Interpreter','latex'); ylabel('$k_x$ (rad/m)','Interpreter','latex');
title('$|S_x(k_x, t)|^2$ ($y$-averaged)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
ylim([0,3000])
ax2 = nexttile;
imagesc(t_plot, kx, log10(Sy_kt + eps));
set(gca, 'YDir', 'normal'); axis tight;
 caxis([0,4])
colormap(ax2, bone);
cb = colorbar; cb.Label.String = '$\log_{10}|S_y(k_x)|^2$';cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;
ylim([0,3000])
xlabel('$t$ (s)','Interpreter','latex'); ylabel('$k_x$ (rad/m)','Interpreter','latex');
title('$|S_y(k_x, t)|^2$ ($y$-averaged)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');

% save_figure(fig, 'kx_t_slopes_csig');

%% Sanity-check figure: per-row spectra before the y-average
j_check = 1;
[~, ~, Px_rows, Py_rows] = row_avg_spectrum( ...
    double(Sx(:,:,idx(j_check))), double(Sy(:,:,idx(j_check))), win_x, Nk);

n_show = min(8, Ny);
y_show = round(linspace(1, Ny, n_show));

fig2 = figure('Name', 'Per-row sanity check', ...
              'Position', [100 100 1200 500], 'Color', 'w');

subplot(1,2,1);
loglog(kx(2:end), Px_rows(y_show, 2:end).', 'Color', [.6 .6 .6]); hold on;
loglog(kx(2:end), mean(Px_rows(:, 2:end), 1), 'k', 'LineWidth', 2);
xlabel('k_x (rad/m)'); ylabel('|S_x(k_x)|^2');
title(sprintf('S_x per-row (grey) vs y-average (black) at t = %.2f s', t_plot(j_check)));
grid on;

subplot(1,2,2);
loglog(kx(2:end), Py_rows(y_show, 2:end).', 'Color', [.6 .6 .6]); hold on;
loglog(kx(2:end), mean(Py_rows(:, 2:end), 1), 'k', 'LineWidth', 2);
xlabel('k_x (rad/m)'); ylabel('|S_y(k_x)|^2');
title(sprintf('S_y per-row (grey) vs y-average (black) at t = %.2f s', t_plot(j_check)));
grid on;
end


function [px_mean, py_mean, Px_rows, Py_rows] = row_avg_spectrum(fx, fy, win_x, Nk)
% Per-row Hann-windowed FFT in x, return y-averaged one-sided power and
% (optionally) the per-row power matrices for sanity plots.

fx_d = fx - mean(fx, 2);                 % remove row mean
fy_d = fy - mean(fy, 2);

Fx = fft(fx_d .* win_x, [], 2);          % FFT along x (dim 2)
Fy = fft(fy_d .* win_x, [], 2);

Px_rows = abs(Fx(:, 1:Nk)).^2;
Py_rows = abs(Fy(:, 1:Nk)).^2;

px_mean = mean(Px_rows, 1).';            % column [Nk x 1]
py_mean = mean(Py_rows, 1).';
end
