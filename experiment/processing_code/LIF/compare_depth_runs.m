%% Compare dye deepening between runs (e.g. stratified 118 vs unstratified 96)
% Each run must have been processed by surf_theta_video.m (run_num = N),
% which writes D:\HLAB_2026\Core\CoreView_N\CoreView_N_surf_theta_z99.mat
% with depths already in cm on that run's own plate calibration
% (util/lif_run_config.m), so runs with different calibrations compare
% directly.  Series are conditioned by util/lif_depth_series.m (1 s running
% median, side-view floor) exactly as in plot_depth_timeseries.m.
%
%   left : depth below the surface vs time since wind start, one colour per
%          run; solid = resolved (deeper than the run's side-view floor),
%          thin grey = within the floor; dotted = each run's floor
%   right: CT profile of the stratified run(s) on the same depth axis
%          (conductivity voltage; the pycnocline -- strongest gradient -- is
%          shaded on both panels: the density gradient refracts light there,
%          so the side-view depth is least reliable in that band)
%   below: hot-wire U of each run, then wave-gauge eta (own gauge or the
%          registry stand-in, e.g. wave_gauge_20261008_231731_r1 for strat)
clearvars -except group run_queue; close all;   % set group before calling (see settings)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'util'));

%% ---------------- User settings ----------------
% Which runs (from processing_code/common/run_registry.csv; only runs with a
% side view; unprocessed ones are skipped with a message):
%   'unstrat' - all unstratified runs  (run-to-run variation)
%   'strat'   - all stratified runs
%   'best'    - the best run of each condition (strat vs unstrat)
%   or a list, e.g. [95 118]
if ~exist('group', 'var'), group = 'best'; end
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'common'));
if ischar(group)
    switch group
        case {'unstrat', 'strat'}, runs = run_registry('list', 'condition', group, 'has', 'side_cam');
        case 'best',               runs = run_registry('list', 'best', 1, 'has', 'side_cam');
    end
else
    runs = group;
end
measure   = 'z99_all';     % 'z99_all' | 'z90_all' | 'zrms_all'
show_rms  = false;         % also draw z_rms (dashed) for each run (busy with many runs)
med_win_s = 1;  jump_cm = 1;  floor_win = [20 40];  floor_fac = 1.25;
pyc_frac  = 0.5;           % pycnocline = |dV/dz| above this fraction of its max
hw_avg_s  = 0.1;
zlim_cm   = [-22 0];
if ischar(group), gname = group; else, gname = strjoin(string(group), '_'); end
out_png   = fullfile('D:\HLAB_2026\Core', ['compare_depth_' char(gname) '.png']);   % '' = don't save

%% ---------------- Which runs are ready ----------------
ok = false(size(runs));
for k = 1:numel(runs)
    ts = fullfile('D:\HLAB_2026\Core', sprintf('CoreView_%d', runs(k)), ...
                  sprintf('CoreView_%d_surf_theta_z99.mat', runs(k)));
    try
        lif_run_config(runs(k));
        ok(k) = isfile(ts);
        if ~ok(k), fprintf('skip run %d: not processed yet (run surf_theta_video / process_runs)\n', runs(k)); end
    catch ME
        fprintf('skip run %d: %s\n', runs(k), ME.message);
    end
end
runs = runs(ok);
if isempty(runs), error('None of the requested runs is processed.'); end

%% ---------------- Load ----------------
nr = numel(runs);
cols = lines(max(nr, 2));
ms = unique([{measure}, {'zrms_all'}], 'stable');
lopt = struct('med_win_s', med_win_s, 'jump_cm', jump_cm, 'floor_win', floor_win, ...
              'floor_fac', floor_fac, 'measures', {ms});
Dr = cell(nr, 1);  Hw = cell(nr, 1);  Wv = cell(nr, 1);  Rr = cell(nr, 1);  Pr = cell(nr, 1);
for k = 1:nr
    Rr{k} = lif_run_config(runs(k));
    ts = fullfile('D:\HLAB_2026\Core', sprintf('CoreView_%d', runs(k)), ...
                  sprintf('CoreView_%d_surf_theta_z99.mat', runs(k)));
    if ~isfile(ts)
        error('No time series for run %d (%s): run surf_theta_video with run_num = %d first.', ...
              runs(k), ts, runs(k));
    end
    Dr{k} = lif_depth_series(ts, lopt);
    [~, Hw{k}, Wv{k}] = lif_wind_clock(Rr{k}, hw_avg_s, fullfile(here, 'calibration'));   % [] = no wind / wave record
    if isempty(Wv{k}) && ~isempty(Hw{k}) && isfield(Hw{k}, 'kind') && strcmp(Hw{k}.kind, 'eta'), Wv{k} = Hw{k}; Hw{k} = []; end   % gauge-only runs (102): eta goes in the wave panel
    M = Dr{k}.(measure);
    fprintf('run %-24s %d frames, floor %.2f cm, resolved from t = %.1f s, final %s = %.2f cm\n', ...
            Rr{k}.label, numel(Dr{k}.t), M.floor, M.t_on, measure, M.filt(find(~isnan(M.filt), 1, 'last')));
    if ~isempty(Rr{k}.density_file)                    % CT profile: depth cm positive down
        Pr{k} = lif_ct_profile(Rr{k}.density_file, pyc_frac);
        fprintf('   CT %s: %d points, pycnocline %.1f..%.1f cm\n', Rr{k}.density_file, numel(Pr{k}.z), Pr{k}.pyc);
    end
end

%% ---------------- Plot ----------------
fig = figure('Color', 'w', 'Position', [80 80 1250 950]);
tl  = tiledlayout(fig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1, [2 3]);  hold(ax1, 'on');
for k = 1:nr                                           % pycnocline bands
    if isempty(Pr{k}), continue; end
    patch(ax1, [0 1e3 1e3 0], Pr{k}.pyc([1 1 2 2]), cols(k,:), 'FaceAlpha', 0.12, ...
          'EdgeColor', 'none', 'DisplayName', sprintf('pycnocline, run %d', runs(k)));
end
for k = 1:nr
    t = Dr{k}.t;  M = Dr{k}.(measure);
    zr = M.filt;  zr(~M.resolved) = NaN;
    plot(ax1, t, M.filt, '-', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
    plot(ax1, t, zr, '-', 'Color', cols(k,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('%s, run %s%s', strrep(measure, '_all', ''), Rr{k}.label, hstr(Rr{k})));
    yline(ax1, M.floor, ':', 'Color', cols(k,:), 'HandleVisibility', 'off');
    if any(Dr{k}.deep)
        plot(ax1, t(Dr{k}.deep), zlim_cm(1) * ones(nnz(Dr{k}.deep), 1), 'v', 'Color', cols(k,:), ...
             'MarkerFaceColor', cols(k,:), 'DisplayName', sprintf('run %d: dye below view', runs(k)));
    end
    if show_rms && ~strcmp(measure, 'zrms_all')
        Q = Dr{k}.zrms_all;  qr = Q.filt;  qr(~Q.resolved) = NaN;
        plot(ax1, t, qr, '--', 'Color', cols(k,:), 'LineWidth', 1.2, ...
             'DisplayName', sprintf('z_{rms}, run %d', runs(k)));
    end
end
hold(ax1, 'off'); grid(ax1, 'on'); box(ax1, 'on');
ylim(ax1, zlim_cm);
ylabel(ax1, 'depth below the surface (cm)');
legend(ax1, 'Location', 'southwest');
title(ax1, sprintf('Dye deepening: %s, %g s running median; grey = within side-view floor', ...
      strrep(measure, '_all', ''), med_win_s), 'Interpreter', 'none');
set(ax1, 'FontSize', 12, 'FontName', 'times');

ax2 = nexttile(tl, 4, [2 1]);  hold(ax2, 'on');       % CT profile(s)
has_ct = false;
for k = 1:nr
    if isempty(Pr{k}), continue; end
    has_ct = true;
    patch(ax2, [min(Pr{k}.v) max(Pr{k}.v) max(Pr{k}.v) min(Pr{k}.v)], Pr{k}.pyc([1 1 2 2]), ...
          cols(k,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax2, Pr{k}.v, Pr{k}.z, 'o-', 'Color', cols(k,:), 'MarkerFaceColor', cols(k,:), ...
         'DisplayName', sprintf('run %d', runs(k)));
end
hold(ax2, 'off'); grid(ax2, 'on'); box(ax2, 'on');
ylim(ax2, zlim_cm);
xlabel(ax2, 'CT conductivity (V)');  title(ax2, 'stratification');
if has_ct, legend(ax2, 'Location', 'southeast'); else, text(ax2, 0.5, 0.5, 'no CT profile', 'Units', 'normalized', 'HorizontalAlignment', 'center'); end
set(ax2, 'FontSize', 12, 'FontName', 'times');
linkaxes([ax1 ax2], 'y');

% Wind (hot-wire / reference probe U) and, below, wave-gauge eta of each run
% (its own gauge, or the stand-in record, labelled as in the videos)
ax3 = nexttile(tl, 9, [1 3]);  hold(ax3, 'on');
for k = find(~cellfun(@isempty, Hw(:)'))
    plot(ax3, Hw{k}.t, Hw{k}.Y(:, 1), '-', 'Color', cols(k,:), 'LineWidth', 0.8, ...
         'DisplayName', sprintf('U, run %d%s', runs(k), src(Hw{k})));
end
hold(ax3, 'off'); grid(ax3, 'on'); box(ax3, 'on');
ylabel(ax3, 'U (m/s)');  set(ax3, 'XTickLabel', []);
legend(ax3, 'Location', 'northwest');
set(ax3, 'FontSize', 12, 'FontName', 'times');

ax4 = nexttile(tl, 13, [1 3]);  hold(ax4, 'on');
for k = find(~cellfun(@isempty, Wv(:)'))
    plot(ax4, Wv{k}.t, Wv{k}.Y(:, 1), '-', 'Color', [cols(k,:) 0.7], 'LineWidth', 0.6, ...
         'DisplayName', sprintf('\\eta, run %d%s', runs(k), src(Wv{k})));
end
hold(ax4, 'off'); grid(ax4, 'on'); box(ax4, 'on');
ylabel(ax4, '\eta (cm)');  xlabel(ax4, 'time since wind start (s)');
if any(~cellfun(@isempty, Wv)), legend(ax4, 'Location', 'northwest');
else, text(ax4, 0.5, 0.5, 'no wave-gauge record', 'Units', 'normalized', 'HorizontalAlignment', 'center'); end
set(ax4, 'FontSize', 12, 'FontName', 'times');
linkaxes([ax1 ax3 ax4], 'x');
tmax = max(cellfun(@(D) D.t(end), Dr));
xlim(ax1, [0 tmax]);

if ~isempty(out_png)
    exportgraphics(fig, out_png, 'Resolution', 200);
    fprintf('Saved %s\n', out_png);
end

%% ---------------- helpers ----------------
function s = src(S)
    % '(simultaneous)' / '(stand-in: core N)' from lif_wind_clock's legend name
    s = regexp(S.names{1}, '\((simultaneous|stand-in:[^)]*)\)', 'match', 'once');
    if ~isempty(s), s = [' ' s]; end                   % file names keep '\_' (tex legend)
end

function s = hstr(R)
    s = '';
    if ~isempty(R.h_cm) && ~contains(R.label, 'h ='), s = sprintf(' (h = %g cm)', R.h_cm); end
end
