%% Transverse vs side-view dye deepening (z99 against time since wind start)
% The paper's LIF plane is transverse (y-z), where the dye is ~uniform
% across the tank; the side view integrates through the tank width.  This
% compares the two for runs on the same day and water:
%   unstrat : side 102  vs transverse 105 (and 104)
%   strat   : side 118  vs transverse 119;  side 120 vs transverse 121
% No run has both views at once, so each pair is two adjacent runs.
%
% Inputs: <run>_surf_theta_z99.mat (side, surf_theta_video.m) and
% <run>_trans_conc_z99.mat (transverse, trans_conc_video.m).  Both are
% conditioned the same way: running median over med_win_s.  Missing runs are
% skipped with a message.
clearvars -except run_queue; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'util'), fullfile(here, '..', 'common'));

%% ---------------- User settings ----------------
pairs     = {[102 105 104], [118 119], [120 121]};   % {side, transverse...}
med_win_s = 1;
zlim_cm   = [-22 0];
core      = 'D:\HLAB_2026\Core';
out_png   = fullfile(core, 'compare_trans_vs_side.png');

side_file  = @(r) fullfile(core, sprintf('CoreView_%d', r), sprintf('CoreView_%d_surf_theta_z99.mat', r));
trans_file = @(r) fullfile(core, sprintf('CoreView_%d', r), sprintf('CoreView_%d_trans_conc_z99.mat', r));

%% ---------------- Plot ----------------
np  = numel(pairs);
fig = figure('Color', 'w', 'Position', [60 60 560*np 520]);
tl  = tiledlayout(fig, 1, np, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:np
    ax = nexttile(tl); hold(ax, 'on');
    runs = pairs{p};  ttl = {};
    % side view
    r = runs(1);
    if isfile(side_file(r))
        D = lif_depth_series(side_file(r), struct('med_win_s', med_win_s, 'measures', {{'z99_all'}}));
        M = D.z99_all;  zr = M.filt;  zr(~M.resolved) = NaN;
        plot(ax, D.t, M.filt, '-', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
        plot(ax, D.t, zr, '-', 'Color', [0 0.45 0.74], 'LineWidth', 2, 'DisplayName', sprintf('side, run %d', r));
        ttl{end+1} = sprintf('side %d', r); %#ok<SAGROW>
    else
        fprintf('skip side run %d: %s not found (run surf_theta_video)\n', r, side_file(r));
    end
    % transverse view(s)
    cols = [0.85 0.33 0.10; 0.93 0.69 0.13];
    for k = 2:numel(runs)
        r = runs(k);
        if ~isfile(trans_file(r))
            fprintf('skip transverse run %d: %s not found (run trans_conc_video)\n', r, trans_file(r));
            continue
        end
        S = load(trans_file(r), 't_wind', 'z99_all', 't_what');
        t = S.t_wind(:);  z = S.z99_all(:);
        kw = max(1, 2*round(med_win_s / median(diff(t)) / 2) + 1);
        lbl = sprintf('transverse, run %d', r);
        if contains(S.t_what, 'assumed'), lbl = [lbl ' (t assumed)']; end
        plot(ax, t, movmedian(z, kw, 'omitnan'), '-', 'Color', cols(k-1,:), 'LineWidth', 2, 'DisplayName', lbl);
        ttl{end+1} = sprintf('trans %d', r); %#ok<SAGROW>
    end
    hold(ax, 'off'); grid(ax, 'on'); box(ax, 'on');
    ylim(ax, zlim_cm);
    xlabel(ax, 'time since wind start (s)');
    if p == 1, ylabel(ax, 'z_{99} below the surface (cm)'); end
    G = run_registry(runs(1));
    title(ax, sprintf('%s: %s', G.condition, strjoin(ttl, ' vs ')));
    legend(ax, 'Location', 'southwest');
    set(ax, 'FontSize', 12, 'FontName', 'times');
end
if ~isempty(out_png)
    exportgraphics(fig, out_png, 'Resolution', 200);
    fprintf('Saved %s\n', out_png);
end
