function plot_wave_ramp_run(src, win_dur, use_omega, clim_log)
% PLOT_WAVE_RAMP_RUN  Re-draw the run_wavegauge_fan_calibration.m summary
% figure (fan command, eta, eta spectrogram, recorded trigger) from a saved
% wave_ramp_*.mat, in the same style as
% Slope_Gauge/plotting/plot_kx_t_slopes_csig.m (tiledlayout, inferno, LaTeX
% labels, Times 14 pt).
%
%   plot_wave_ramp_run('wave_ramp_20260928_120000.mat')
%   plot_wave_ramp_run(data)                    % struct already loaded with load()
%   plot_wave_ramp_run(file, 4, false)          % f (Hz) axis instead of omega
%   plot_wave_ramp_run(file, 4, true, [-8 -4])  % fixed colour limits (log10)
%
% win_dur   : spectrogram window (s), default 4
% use_omega : true (default) = omega in rad/s, S_eta(omega) = S_eta(f)/(2 pi),
%             consistent with k in rad/m for the later k-omega combination
%             with the CSIG slopes; false = f in Hz
% clim_log  : colour limits for log10 S_eta ([] = auto; the CSIG plots fix
%             theirs with caxis([0 4]))

if nargin < 2 || isempty(win_dur),   win_dur = 4;      end
if nargin < 3 || isempty(use_omega), use_omega = true; end
if nargin < 4, clim_log = []; end

% inferno lives with the CSIG code
SLOPE_UTIL = 'D:\Chris\osbl-turbulent-mixing\experiment\processing_code\Slope_Gauge\util';
if ~exist('inferno', 'file') && exist(SLOPE_UTIL, 'dir'), addpath(SLOPE_UTIL); end

if ischar(src) || isstring(src)
    [~, fname] = fileparts(src);
    d = load(src);
else
    d = src; fname = 'wave ramp run';
end

FS      = d.runConfig.FS;
logTrig = isfield(d, 'E_trig') && ~isempty(d.E_trig);
camOnT  = d.camOnT;  camOffT = d.camOffT;
F_MAX   = 10;   % Hz, display limit (wire gauge resolution)
etaU    = d.etaUnit;

nSub = 3 + double(logTrig);
figure('Name', sprintf('Wave ramp run: %s', fname), ...
       'Position', [100 100 1100 250*nSub], 'Color', 'w');
tiledlayout(nSub, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(nSub, 1);

%% Fan command
ax(1) = nexttile;
% plot (not stairs), as plot_experiment_signals.m draws run_experiment.m's
% fan -- the AO holds each write for FAN_DT, which stairs would show as
% ~0.9 V steps on the 5 s ramp-down.
plot(d.sentT - d.fanStartElapsed, d.sentV, 'LineWidth', 1.2);
shade_camera(camOnT, camOffT);
ylabel('Fan (V)', 'Interpreter', 'latex');
title(sprintf('Fan command, %s (shaded = camera triggering)', strrep(fname, '_', '\_')), ...
      'Interpreter', 'latex');
style_axes(ax(1), true);

%% eta(t)
ax(2) = nexttile;
plot(d.hwT_wind, d.eta, 'LineWidth', 0.5);
shade_camera(camOnT, camOffT);
ylabel(sprintf('$\\eta$ (%s)', etaU), 'Interpreter', 'latex');
title('Wave gauge', 'Interpreter', 'latex');
style_axes(ax(2), true);

%% Spectrogram
ax(3) = nexttile;
% 4th output of spectrogram is the one-sided PSD in etaUnit^2/Hz.
nWin = round(win_dur * FS);
[~, fS, tS, pS] = spectrogram(detrend(d.eta), hann(nWin), round(nWin/2), [], FS);
if use_omega
    yS = 2*pi * fS;  pS = pS / (2*pi);   % S(omega) = S(f)/(2 pi): same variance
    yLab  = '$\omega$ (rad/s)';  yMax = 2*pi * F_MAX;
    cbLab = sprintf('$\\log_{10} S_\\eta(\\omega)$ (%s$^2$ s/rad)', etaU);
    ttl   = sprintf('$S_\\eta(\\omega, t)$ (%g s windows)', win_dur);
else
    yS = fS;
    yLab  = '$f$ (Hz)';  yMax = F_MAX;
    cbLab = sprintf('$\\log_{10} S_\\eta(f)$ (%s$^2$/Hz)', etaU);
    ttl   = sprintf('$S_\\eta(f, t)$ (%g s windows)', win_dur);
end
imagesc(tS + d.hwT_wind(1), yS, log10(pS + eps));
set(gca, 'YDir', 'normal'); axis tight;
if ~isempty(clim_log), caxis(clim_log); end
colormap(ax(3), pick_cmap());
cb = colorbar;
cb.Label.String = cbLab; cb.Label.Interpreter = 'latex'; cb.Label.FontSize = 16;
ylim([0 yMax]);
ylabel(yLab, 'Interpreter', 'latex');
title(ttl, 'Interpreter', 'latex');
shade_camera(camOnT, camOffT);
style_axes(ax(3), false);

%% Recorded trigger
if logTrig
    ax(4) = nexttile;
    plot(d.hwT_wind, d.E_trig); hold on;
    plot(d.frameT, d.runConfig.CAM_TRIG_THRESH * ones(size(d.frameT)), 'r.', 'MarkerSize', 4); hold off;
    ylabel('Trigger (V)', 'Interpreter', 'latex');
    title('Recorded camera trigger (red = detected frames)', 'Interpreter', 'latex');
    style_axes(ax(4), true);
end
xlabel('$t$ (s)', 'Interpreter', 'latex');

% Invisible colorbars on the line-plot tiles so every axes has the same
% width as the spectrogram and the linked time axes line up.
for i = [1 2 4]
    if i <= nSub
        c = colorbar(ax(i)); c.Visible = 'off';
    end
end
linkaxes(ax, 'x');
xlim(ax(1), [d.hwT_wind(1), d.hwT_wind(end)]);
end


function style_axes(a, withGrid)
% Same font settings as the CSIG plots.
    set(a, 'fontsize', 14, 'fontname', 'times');
    if withGrid, grid(a, 'on'); end
end

function cm = pick_cmap()
    if exist('inferno', 'file'), cm = inferno; else, cm = parula; end
end

function shade_camera(tOn, tOff)
% Magenta band behind each camera window (same as run_wavegauge_fan_calibration.m).
    if isempty(tOn), return; end
    wasHeld = ishold; hold on;
    try
        h = xregion(tOn(:), tOff(:), 'FaceColor', [0.85 0.2 0.85], 'FaceAlpha', 0.12);
        set(h, 'HandleVisibility', 'off');
    catch
        yl = ylim;
        for k = 1:numel(tOn)
            patch([tOn(k) tOff(k) tOff(k) tOn(k)], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.85 0.2 0.85], 'EdgeColor','none', 'FaceAlpha',0.12, 'HandleVisibility','off');
        end
        ylim(yl);
    end
    if ~wasHeld, hold off; end
end
