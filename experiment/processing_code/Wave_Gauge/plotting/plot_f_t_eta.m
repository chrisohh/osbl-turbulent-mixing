function S = plot_f_t_eta(data, win_dur, center_times, t_range)
% PLOT_F_T_ETA  Frequency-vs-time spectrogram of the wave-gauge elevation.
% Gauge counterpart of Slope_Gauge/plotting/plot_f_t_slopes_csig.m: same
% figures (f-t map with peak ridge, peak f / peak power vs t, spectra at
% chosen times), but computed directly from eta(t) -- no dispersion mapping
% is needed because the gauge already measures in time.
%
%   data = load('wave_ramp_20260928_120000.mat');
%   plot_f_t_eta(data)                      % 4 s windows, whole record
%   plot_f_t_eta(data, 4, [20 40 60], 'cam')
%
% win_dur      : window length (s), default 4
% center_times : times (s, wind clock) for the per-time spectra figure
%                (default: 6 evenly spaced)
% t_range      : [t1 t2] or 'cam' (camera window only), default whole record

if nargin < 2 || isempty(win_dur), win_dur = 4; end
if nargin < 3, center_times = []; end
if nargin < 4, t_range = []; end

S = eta_spectrogram(data, win_dur, t_range);

F_LO = 0.5;   % Hz, ignore drift / seiching when picking the peak (as the slope code)
iLo  = find(S.f >= F_LO, 1);
[S.PPeak, iPk] = max(S.P(iLo:end, :), [], 1);
S.fPeak = S.f(iPk + iLo - 1).';

cbLabel = sprintf('$\\log_{10} S_\\eta(f)$ (%s$^2$/Hz)', S.unit);

%% f-t spectrogram
figure('Name', 'f-t spectrogram (eta)', 'Position', [100 100 1100 450], 'Color', 'w');
imagesc(S.t, S.f, log10(S.P + eps));
set(gca, 'YDir', 'normal'); axis tight;
colormap(gca, pick_cmap());
cb = colorbar;
cb.Label.String = cbLabel; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;
hold on;
plot(S.t, S.fPeak, 'r-', 'LineWidth', 1.5);
mark_camera(S);
xlabel('$t$ (s)',  'Interpreter', 'latex');
ylabel('$f$ (Hz)', 'Interpreter', 'latex');
title(sprintf('$S_\\eta(f,t)$, wave gauge (%g s windows)', S.win_dur), 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times');

%% Peak frequency and peak power evolution
figure('Name', 'Peak frequency and power (eta)', 'Position', [100 100 1000 600], 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(S.t, S.fPeak, 'b-', 'LineWidth', 1.5); hold on; mark_camera(S);
ylabel('$f_p$ (Hz)', 'Interpreter', 'latex');
title('Peak frequency', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times'); grid on;

nexttile;
semilogy(S.t, S.PPeak, 'b-', 'LineWidth', 1.5); hold on; mark_camera(S);
xlabel('$t$ (s)', 'Interpreter', 'latex');
ylabel(sprintf('$S_\\eta(f_p)$ (%s$^2$/Hz)', S.unit), 'Interpreter', 'latex');
title('Spectral power at the peak', 'Interpreter', 'latex');
set(gca, 'fontsize', 14, 'fontname', 'times'); grid on;

%% Spectra at selected times
if isempty(center_times)
    center_times = linspace(S.t(1), S.t(end), 6);
end
[~, jc] = min(abs(S.t(:) - center_times(:).'), [], 1);
colors  = lines(numel(jc));
leg_str = arrayfun(@(t) sprintf('t = %.0f s', t), S.t(jc), 'UniformOutput', false);

figure('Name', 'eta spectra at each time', 'Position', [100 100 700 500], 'Color', 'w');
for j = 1:numel(jc)
    loglog(S.f(2:end), S.P(2:end, jc(j)), 'Color', colors(j,:), 'LineWidth', 1.2); hold on;
end
xlabel('$f$ (Hz)', 'Interpreter', 'latex');
ylabel(sprintf('$S_\\eta(f)$ (%s$^2$/Hz)', S.unit), 'Interpreter', 'latex');
title('$\eta$ spectra', 'Interpreter', 'latex');
legend(leg_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 9);
set(gca, 'fontsize', 12, 'fontname', 'times'); grid on;
end


function cm = pick_cmap()
% inferno as in the CSIG plots when it is on the path, parula otherwise.
    if exist('inferno', 'file'), cm = inferno; else, cm = parula; end
end

function mark_camera(S)
    if ~isnan(S.camOn)
        xline([S.camOn S.camOff], 'm--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
