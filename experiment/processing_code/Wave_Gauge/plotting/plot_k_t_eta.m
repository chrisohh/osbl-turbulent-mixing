function K = plot_k_t_eta(data, win_dur, center_times, t_range)
% PLOT_K_T_ETA  Wavenumber-vs-time spectrogram from the wave gauge.
% Gauge counterpart of Slope_Gauge/plotting/plot_kx_t_slopes_csig.m. The
% gauge is a single point, so k is NOT measured: the f-t spectrum of eta is
% mapped to k through the same linear gravity-capillary dispersion relation
% plot_f_t_slopes_csig.m uses (in reverse),
%     omega = sqrt(g*k + gamma*k^3),   S(k) = S(f) * df/dk,
% then interpolated onto a uniform k grid.
%
% Two panels:
%   top    S_eta(k,t)          elevation spectrum
%   bottom k^2 S_eta(k,t)      slope spectrum the gauge PREDICTS -- same
%                              quantity as the CSIG Sx(kx,t) panel, for
%                              comparing peak k and spectral shape. The CSIG
%                              plots show un-normalised |FFT|^2, so compare
%                              shapes and peak positions, not absolute level.
%
% Assumptions behind the mapping (check before over-interpreting):
%   - linear, deep-water waves propagating along x (k here ~ kx)
%   - no current: wind drift Doppler-shifts f, so k is biased where the drift
%     is comparable to the phase speed (short waves)
%   - bound harmonics of steep waves appear at 2f but are mapped as if free
%
%   data = load('wave_ramp_20260928_120000.mat');
%   plot_k_t_eta(data, 4, [], 'cam')

if nargin < 2 || isempty(win_dur), win_dur = 4; end
if nargin < 3, center_times = []; end
if nargin < 4, t_range = []; end

S = eta_spectrogram(data, win_dur, t_range);

g = 9.81; gamma = 7.4e-5;   % same constants as plot_f_t_slopes_csig.m

% Invert f(k) numerically (monotonic), drop f = 0.
kg  = logspace(-2, 4, 20000);
fkg = sqrt(g*kg + gamma*kg.^3) / (2*pi);
use = S.f > 0;
f   = S.f(use);
k_f = interp1(fkg, kg, f);                                   % [Nf x 1]
dfdk = (g + 3*gamma*k_f.^2) ./ (4*pi*sqrt(g*k_f + gamma*k_f.^3));
Pk_nonuni = S.P(use, :) .* dfdk;                             % per rad/m

% Uniform k grid over the same range (as the f-t slope code does for f)
K.k   = linspace(k_f(1), k_f(end), numel(k_f)).';
K.Pk  = interp1(k_f, Pk_nonuni, K.k, 'linear', 0);           % S_eta(k,t)
K.Bk  = (K.k.^2) .* K.Pk;                                    % k^2 S_eta(k,t)
K.t   = S.t;

slopeOK = strcmp(S.unit, 'm');
if ~slopeOK
    warning('Gauge not calibrated (eta in %s) -- k^2 S panel is not a slope spectrum.', S.unit);
end

%% k-t spectrogram
figure('Name', 'k-t spectrogram (eta, via dispersion)', 'Position', [100 100 1100 800], 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
imagesc(K.t, K.k, log10(K.Pk + eps));
set(gca, 'YDir', 'normal'); axis tight; colormap(ax1, pick_cmap());
cb = colorbar;
cb.Label.String = sprintf('$\\log_{10} S_\\eta(k)$ (%s$^2$ m)', S.unit);
cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;
hold on; mark_camera(S);
ylabel('$k$ (rad/m)', 'Interpreter', 'latex');
title('$S_\eta(k,t)$, wave gauge mapped through dispersion', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');

ax2 = nexttile;
imagesc(K.t, K.k, log10(K.Bk + eps));
set(gca, 'YDir', 'normal'); axis tight; colormap(ax2, pick_cmap());
cb = colorbar;
cb.Label.String = '$\log_{10} k^2 S_\eta(k)$ (rad$^2$ m$^{-1}$)';
cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;
hold on; mark_camera(S);
xlabel('$t$ (s)', 'Interpreter', 'latex'); ylabel('$k$ (rad/m)', 'Interpreter', 'latex');
title('$k^2 S_\eta(k,t)$ -- slope spectrum predicted by the gauge', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');
linkaxes([ax1 ax2], 'xy');

%% Spectra at selected times
if isempty(center_times)
    center_times = linspace(K.t(1), K.t(end), 6);
end
[~, jc] = min(abs(K.t(:) - center_times(:).'), [], 1);
colors  = lines(numel(jc));
leg_str = arrayfun(@(t) sprintf('t = %.0f s', t), K.t(jc), 'UniformOutput', false);

figure('Name', 'k spectra at each time (eta)', 'Position', [100 100 1200 500], 'Color', 'w');
subplot(1,2,1);
for j = 1:numel(jc)
    loglog(K.k, K.Pk(:, jc(j)), 'Color', colors(j,:), 'LineWidth', 1.2); hold on;
end
xlabel('$k$ (rad/m)', 'Interpreter', 'latex');
ylabel(sprintf('$S_\\eta(k)$ (%s$^2$ m)', S.unit), 'Interpreter', 'latex');
title('$\eta$ spectra', 'Interpreter', 'latex');
legend(leg_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
set(gca, 'fontsize', 12, 'fontname', 'times'); grid on;

subplot(1,2,2);
for j = 1:numel(jc)
    loglog(K.k, K.Bk(:, jc(j)), 'Color', colors(j,:), 'LineWidth', 1.2); hold on;
end
xlabel('$k$ (rad/m)', 'Interpreter', 'latex');
ylabel('$k^2 S_\eta(k)$', 'Interpreter', 'latex');
title('Gauge-predicted slope spectra', 'Interpreter', 'latex');
legend(leg_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
set(gca, 'fontsize', 12, 'fontname', 'times'); grid on;
end


function cm = pick_cmap()
    if exist('inferno', 'file'), cm = inferno; else, cm = parula; end
end

function mark_camera(S)
    if ~isnan(S.camOn)
        xline([S.camOn S.camOff], 'm--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
