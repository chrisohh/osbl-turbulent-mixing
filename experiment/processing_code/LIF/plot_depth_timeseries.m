%% Dye depth vs time since wind start, with the side-view resolution floor
% Reads the time series saved by surf_theta_video.m (<video>_z99.mat):
% z99 (paper eq. 3.14), z90 and the RMS dye depth, each the median of the
% 5 cm x-bins (util/lif_z99.m).
%
% Resolution floor: the side view looks THROUGH the ~50 cm tank width, so
% absorption is integrated along the line of sight.  Waves across the width
% and the camera not sitting exactly at the waterline smear a thin surface
% layer into an apparent band of ~1-2 cm.  While the dye is still a surface
% layer every measure therefore sits on a plateau (CoreView 95: t ~ 20-40 s,
% z99 ~ -2 cm, z90 ~ -1 cm, z_rms ~ -0.7 cm) that is the floor of the
% method, not a mixed layer.  Each measure's floor is its median over
% floor_win; points not deeper than floor_fac * floor are drawn grey as
% unresolved, and the band above the z99 floor is shaded.
clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'util'));

%% ---------------- User settings ----------------
ts_file   = 'D:\HLAB_2026\Core\CoreView_95\CoreView_95_surf_theta_z99.mat';
hw_file   = 'D:\HLAB_2026\hotwire\hotwire_20260923_113939.mat';
hw_avg_s  = 0.1;           % hot-wire block average for the lower panel (s)

floor_win = [20 40];       % s since wind start: plateau = surface layer only
floor_fac = 1.25;          % resolved = deeper than floor_fac * floor
measures  = {'z99_all', '$z_{99}$';
             'z90_all', '$z_{90}$';
             'zrms_all', '$z_{\rm rms}$'};
cols      = lines(3);

%% ---------------- Load ----------------
S = load(ts_file);
t = S.t_wind(:);
fprintf('%s: %d frames, t = %.1f..%.1f s\n', ts_file, numel(t), t(1), t(end));

[~, hw_name] = fileparts(hw_file);
[~, hw] = lif_hw_cached(hw_file, hw_avg_s, fullfile(here, 'calibration', ...
                        sprintf('hwsync_%s_avg%g.mat', hw_name, hw_avg_s)));

%% ---------------- Floors ----------------
in_win = t >= floor_win(1) & t <= floor_win(2);
nm = size(measures, 1);
flo = nan(nm, 1);
for m = 1:nm
    z = S.(measures{m,1})(:);
    flo(m) = median(z(in_win), 'omitnan');
    res = z < floor_fac * flo(m);
    t_on = t(find(res & t > floor_win(2), 1, 'first'));
    fprintf('%-8s floor %.2f cm (median over %g-%g s); first resolved after the window: t = %.1f s\n', ...
            measures{m,1}, flo(m), floor_win, t_on);
end

%% ---------------- Plot ----------------
fig = figure('Color', 'w', 'Position', [80 80 1100 750]);
tl  = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, [2 1]); hold(ax1, 'on');
yl = [floor(1.1 * min(S.z99_all, [], 'omitnan')) 0];
patch(ax1, t([1 end end 1]), [flo(1) flo(1) 0 0], [0.85 0.85 0.85], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', 'unresolved (side-view smear, $z_{99}$ floor)');
if isfield(S, 'z99_x_all')                     % spread of the local z99(x)
    q = prctile(S.z99_x_all, [5 95], 2);
    ok = all(isfinite(q), 2);
    fill(ax1, [t(ok); flipud(t(ok))], [q(ok,1); flipud(q(ok,2))], cols(1,:), ...
         'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', '$z_{99}(x)$ 5--95\%');
end
for m = 1:nm
    z = S.(measures{m,1})(:);
    res = z < floor_fac * flo(m);
    zr = z;  zr(~res) = NaN;                   % resolved
    plot(ax1, t, z, '-', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
    plot(ax1, t, zr, '.-', 'Color', cols(m,:), 'LineWidth', 1.5, 'MarkerSize', 10, ...
         'DisplayName', measures{m,2});
    yline(ax1, flo(m), ':', 'Color', cols(m,:), 'HandleVisibility', 'off');
end
xline(ax1, floor_win, '--k', 'HandleVisibility', 'off');
hold(ax1, 'off'); grid(ax1, 'on'); box(ax1, 'on');
ylim(ax1, yl);
ylabel(ax1, 'depth below the surface (cm)');
legend(ax1, 'Interpreter', 'latex', 'Location', 'southwest');
title(ax1, sprintf(['Dye depth; grey = within the side-view floor (median over %g--%g s, ' ...
      'dotted), colour = resolved'], floor_win), 'Interpreter', 'latex');
set(ax1, 'FontSize', 12, 'FontName', 'times');

ax2 = nexttile(tl);
plot(ax2, hw.t, hw.Y, 'LineWidth', 0.8);
grid(ax2, 'on'); box(ax2, 'on');
xlabel(ax2, 'time since wind start (s)'); ylabel(ax2, 'velocity (m/s)');
legend(ax2, hw.names, 'Interpreter', 'latex', 'Location', 'northwest');
set(ax2, 'FontSize', 12, 'FontName', 'times');
linkaxes([ax1 ax2], 'x');
xlim(ax1, [min(t(1), hw.t(1)) max(t(end), hw.t(end))]);
