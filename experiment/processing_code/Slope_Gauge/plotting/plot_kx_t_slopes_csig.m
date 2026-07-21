function plot_kx_t_slopes_csig(data, frame_indices, win_dur, center_indices)
% PLOT_KX_T_SLOPES_CSIG  kx-vs-time spectrogram of Sx and Sy.
%
% frame_indices : centers for the spectrogram ([] = all frames)
% win_dur       : averaging window duration in seconds (0 = single frame)
% center_indices: indices into data.time for the per-time spectra line plot.
%                 If omitted, uses frame_indices (or all frames if that is []).

Sx   = data.Sx;
Sy   = data.Sy;
time = data.time;
dx   = data.dx;

[Ny, Nx, Nt] = size(Sx);

if nargin < 2 || isempty(frame_indices)
    idx = 1:Nt;
else
    idx = frame_indices(:).';
end
n_t = numel(idx);

if nargin < 3 || isempty(win_dur) || win_dur <= 0
    win_dur = 0;
end

if nargin < 4 || isempty(center_indices)
    center_indices = idx;
end

% One-sided wavenumber axis (rad/m).
Nk = floor(Nx/2) + 1;
kx = 2*pi * (0:Nk-1) / (Nx*dx);

win_x = hann(Nx).';

Sx_kt = zeros(Nk, n_t);
Sy_kt = zeros(Nk, n_t);

for j = 1:n_t
    if win_dur > 0
        % frames within ±half-window of the center time
        t_c   = time(idx(j));
        in_win = find(abs(time - t_c) <= win_dur/2);
    else
        in_win = idx(j);
    end
    acc_x = zeros(Nk, 1);
    acc_y = zeros(Nk, 1);
    for k = in_win(:).'
        [px, py] = row_avg_spectrum( ...
            double(Sx(:,:,k)), double(Sy(:,:,k)), win_x, Nk);
        acc_x = acc_x + px;
        acc_y = acc_y + py;
    end
    Sx_kt(:, j) = acc_x / numel(in_win);
    Sy_kt(:, j) = acc_y / numel(in_win);
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
colormap(ax1, inferno);
cb = colorbar; cb.Label.String = '$\log_{10}|S_x(k_x)|^2$';cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;

xlabel('$t$ (s)','Interpreter','latex'); ylabel('$k_x$ (rad/m)','Interpreter','latex');
title('$|S_x(k_x, t)|^2$ ($y$-averaged)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
ylim([0,3000])
ax2 = nexttile;
imagesc(t_plot, kx, log10(Sy_kt + eps));
set(gca, 'YDir', 'normal'); axis tight;
 caxis([0,4])
colormap(ax2, inferno);
cb = colorbar; cb.Label.String = '$\log_{10}|S_y(k_x)|^2$';cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;
ylim([0,3000])
xlabel('$t$ (s)','Interpreter','latex'); ylabel('$k_x$ (rad/m)','Interpreter','latex');
title('$|S_y(k_x, t)|^2$ ($y$-averaged)', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');

% save_figure(fig, 'kx_t_slopes_csig');

%% Per-time spectra: compute for center_indices
n_c = numel(center_indices);
Sx_kc = zeros(Nk, n_c);
Sy_kc = zeros(Nk, n_c);
for j = 1:n_c
    if win_dur > 0
        t_c    = time(center_indices(j));
        in_win = find(abs(time - t_c) <= win_dur/2);
    else
        in_win = center_indices(j);
    end
    acc_x = zeros(Nk, 1);
    acc_y = zeros(Nk, 1);
    for k = in_win(:).'
        [px, py] = row_avg_spectrum( ...
            double(Sx(:,:,k)), double(Sy(:,:,k)), win_x, Nk);
        acc_x = acc_x + px;
        acc_y = acc_y + py;
    end
    Sx_kc(:,j) = acc_x / numel(in_win);
    Sy_kc(:,j) = acc_y / numel(in_win);
end
t_centers = time(center_indices);

fig2 = figure('Name', 'kx spectra at each time', ...
              'Position', [100 100 1200 500], 'Color', 'w');

colors  = lines(n_c);
leg_str = arrayfun(@(t) sprintf('t = %.0f s', t), t_centers, 'UniformOutput', false);

subplot(1,2,1); hold on;
for j = 1:n_c
    loglog(kx(2:end), Sx_kc(2:end, j), 'Color', colors(j,:), 'LineWidth', 1.2);
end
xlabel('$k_x$ (rad/m)', 'Interpreter', 'latex');
ylabel('$|S_x(k_x)|^2$', 'Interpreter', 'latex');
title('$S_x$ spectra ($y$-averaged)', 'Interpreter', 'latex');
legend(leg_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
set(gca, 'fontsize', 12, 'fontname', 'times');
grid on;

subplot(1,2,2); hold on;
for j = 1:n_c
    loglog(kx(2:end), Sy_kc(2:end, j), 'Color', colors(j,:), 'LineWidth', 1.2);
end
xlabel('$k_x$ (rad/m)', 'Interpreter', 'latex');
ylabel('$|S_y(k_x)|^2$', 'Interpreter', 'latex');
title('$S_y$ spectra ($y$-averaged)', 'Interpreter', 'latex');
legend(leg_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
set(gca, 'fontsize', 12, 'fontname', 'times');
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
