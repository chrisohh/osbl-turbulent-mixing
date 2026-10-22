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
clearvars -except run_num; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'util'));

%% ---------------- User settings ----------------
if ~exist('run_num', 'var'), run_num = 95; end   % util/lif_run_config.m
R = lif_run_config(run_num);
ts_file  = fullfile('D:\HLAB_2026\Core', sprintf('CoreView_%d', run_num), sprintf('CoreView_%d_surf_theta_z99.mat', run_num));
hw_file  = R.hw_file;
hw_avg_s  = 0.1;           % hot-wire block average for the lower panel (s)

floor_win = [20 40];       % s since wind start: plateau = surface layer only
floor_fac = 1.25;          % resolved = deeper than floor_fac * floor
measures  = {'z99_all', '$z_{99}$';
             'z90_all', '$z_{90}$';
             'zrms_all', '$z_{\rm rms}$'};
cols      = lines(3);

% Time filtering.  Raw per-frame depths are kept and plotted; the line is a
% running median over med_win_s.  Frames further than jump_cm from it are
% flagged: the dye deepens at most ~0.5 cm/s here, so a jump of ~1 cm
% between frames tenths of a second apart is not deepening (moving faint
% leftover dye, a bad frame).
med_win_s = 1;             % s
jump_cm   = 1;             % flag |raw - running median| > jump_cm

%% ---------------- Load + condition (util/lif_depth_series) ----------------
D = lif_depth_series(ts_file, struct('med_win_s', med_win_s, 'jump_cm', jump_cm, ...
        'floor_win', floor_win, 'floor_fac', floor_fac, 'measures', {measures(:,1)'}));
S = load(ts_file);
t = D.t;
fprintf('%s: %d frames, t = %.1f..%.1f s\n', ts_file, numel(t), t(1), t(end));

[~, hw_name] = fileparts(hw_file);
[~, hw] = lif_hw_cached(hw_file, hw_avg_s, fullfile(here, 'calibration', ...
                        sprintf('hwsync_%s_avg%g.mat', hw_name, hw_avg_s)));

nm = size(measures, 1);
Zraw = cell(nm, 1);  Zf = cell(nm, 1);  flag = cell(nm, 1);  flo = nan(nm, 1);
for m = 1:nm
    M = D.(measures{m,1});
    Zraw{m} = M.raw;  Zf{m} = M.filt;  flag{m} = M.flag;  flo(m) = M.floor;
    fprintf('%-8s running median over %d frames; %d of %d frames flagged (> %g cm off); floor %.2f cm, resolved from t = %.1f s\n', ...
            measures{m,1}, D.kw, nnz(M.flag), nnz(~isnan(M.raw)), jump_cm, M.floor, M.t_on);
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
    z = Zf{m};                                 % running median
    res = z < floor_fac * flo(m);
    zr = z;  zr(~res) = NaN;                   % resolved
    plot(ax1, t, Zraw{m}, '.', 'Color', 0.5*cols(m,:) + 0.5, 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(ax1, t(flag{m}), Zraw{m}(flag{m}), 'x', 'Color', [0.85 0 0], 'MarkerSize', 6, ...
         'HandleVisibility', 'off');
    plot(ax1, t, z, '-', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
    plot(ax1, t, zr, '-', 'Color', cols(m,:), 'LineWidth', 1.8, ...
         'DisplayName', [measures{m,2} ' (running median)']);
    yline(ax1, flo(m), ':', 'Color', cols(m,:), 'HandleVisibility', 'off');
end
xline(ax1, floor_win, '--k', 'HandleVisibility', 'off');
% Late frames: dye reaches the bottom of the view in most of x, the edge is
% not visible and z99 is undefined (NaN) -- shade those times; z90 / z_rms
% there are lower bounds
if isfield(S, 'deep_all')
    dp = S.deep_all(:) > 0.5;
    if any(dp)
        e = diff([0; dp; 0]);  s0 = find(e == 1);  s1 = find(e == -1) - 1;
        for k = 1:numel(s0)
            patch(ax1, t([s0(k) s1(k) s1(k) s0(k)]), yl([1 1 2 2]), [1 0.85 0.7], ...
                  'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
        end
        plot(ax1, NaN, NaN, 's', 'MarkerFaceColor', [1 0.85 0.7], 'MarkerEdgeColor', 'none', ...
             'DisplayName', 'dye below the field of view ($z_{99}$ undefined)');
        fprintf('dye below the field of view in %d frames, from t = %.1f s\n', nnz(dp), t(find(dp, 1)));
    end
end
hold(ax1, 'off'); grid(ax1, 'on'); box(ax1, 'on');
ylim(ax1, yl);
ylabel(ax1, 'depth below the surface (cm)');
legend(ax1, 'Interpreter', 'latex', 'Location', 'southwest');
title(ax1, sprintf(['Dye depth: dots = per frame, line = %g s running median, ' ...
      'red x = jump $>$ %g cm; grey = within the side-view floor'], med_win_s, jump_cm), ...
      'Interpreter', 'latex');
set(ax1, 'FontSize', 12, 'FontName', 'times');

ax2 = nexttile(tl);
plot(ax2, hw.t, hw.Y, 'LineWidth', 0.8);
grid(ax2, 'on'); box(ax2, 'on');
xlabel(ax2, 'time since wind start (s)'); ylabel(ax2, 'velocity (m/s)');
legend(ax2, hw.names, 'Interpreter', 'latex', 'Location', 'northwest');
set(ax2, 'FontSize', 12, 'FontName', 'times');
linkaxes([ax1 ax2], 'x');
xlim(ax1, [min(t(1), hw.t(1)) max(t(end), hw.t(end))]);
