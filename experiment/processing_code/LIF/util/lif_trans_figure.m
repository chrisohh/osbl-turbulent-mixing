function h = lif_trans_figure(y_cm, z_cm, opts)
%LIF_TRANS_FIGURE Three-panel transverse-LIF figure: raw | C | profile.
%
%   h = LIF_TRANS_FIGURE(y_cm, z_cm, opts)   builds the (empty) figure once
%   LIF_TRANS_UPDATE(h, ...)                 fills it for one frame
%
% The transverse (y, z) counterpart of lif_surf_figure, shared by
% trans_frame_calibrated.m (one frame) and trans_conc_video.m (every
% frame) so the still and the video look the same.  Differences from the
% side-view figure, all because the dye emits rather than absorbs:
%   * the colour map is NOT flipped -- dye is bright, like the raw image;
%   * the raw panel is log-scaled by default, since the surface line is
%     ~150 counts and the deep plumes ~1 count over a black background;
%   * the overlays are the surface and z99(y), plus the dye-penetration
%     contour z_pen(y) from lif_trans_conc.
%
%   1. raw counts (opts.raw_log: log colour scale)
%   2. concentration C from lif_trans_conc, in the units set by its
%      normalize option
%   3. y-averaged C vs depth below the local surface, z99 and the 5-95%
%      spread of z99(y)
%
% Layout: the transverse field is ~25 cm wide and ~5 cm deep, so at a true
% 1:1 aspect the two images are long letterbox strips -- they are STACKED
% (raw above C) with the profile down the right-hand side, rather than set
% side by side as in the near-square side view.
%
% opts fields (defaults):
%   cmap 'inferno' (falls back to 'hot' if not on the path), clims [0 150],
%   raw_log true, clims_log [0.5 300], C_clims [0 1], C_log false,
%   C_clims_log [1e-3 1], C_label '$C$ (normalised)',
%   profile_xlim [], profile_zlim [-4 0.5], fig_pos [40 40 2000 650],
%   visible 'on',
%   zlim_cm [] / ylim_cm [] -- image-panel limits ([] = whole frame).  The
%      analysis wants the WHOLE frame (lif_trans_conc needs clear water
%      well below the dye for its noise floor), while the dye occupies the
%      top fifth of it, so zoom here rather than cropping the data,
%   hw [] -- hot-wire panel under the images (full width) with a vertical
%      line at the frame time, exactly as in lif_surf_figure.  Struct with
%      t (s since wind start), Y ([nt x nseries] m/s), names (cell, LaTeX),
%      xlim ([] = whole record).

    def = struct('cmap', 'inferno', 'clims', [0 150], ...
                 'raw_log', true, 'clims_log', [0.5 300], ...
                 'C_clims', [0 1], 'C_log', false, 'C_clims_log', [1e-3 1], ...
                 'C_label', '$C$ (normalised)', ...
                 'profile_xlim', [], 'profile_zlim', [-4 0.5], ...
                 'fig_pos', [40 40 1700 760], 'visible', 'on', ...
                 'zlim_cm', [], 'ylim_cm', [], 'hw', []);
    if nargin < 3, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end
    if ischar(opts.cmap) && isempty(which(opts.cmap)) && ~any(strcmp(opts.cmap, ...
            {'parula','jet','hot','bone','gray','hsv','cool','spring','summer','autumn','winter','copper','pink','turbo'}))
        opts.cmap = 'hot';
    end
    cmapM = feval(opts.cmap, 256);
    c99 = [0.3 0.75 1];  cpen = [0.4 1 0.5];
    nx = numel(y_cm);  ny = numel(z_cm);

    has_hw = ~isempty(opts.hw);
    if has_hw && isequal(opts.fig_pos, def.fig_pos), opts.fig_pos(4) = 1050; end
    h.fig = figure('Position', opts.fig_pos, 'Color', 'w', 'Visible', opts.visible);
    % 6 columns: images span 5, the profile the 6th over both image rows.
    tl = tiledlayout(h.fig, 2 + has_hw, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 1. raw counts
    h.axR = nexttile(tl, 1, [1 5]);
    h.imR = imagesc(h.axR, y_cm, z_cm, nan(ny, nx));
    colormap(h.axR, cmapM);
    if opts.raw_log
        set(h.axR, 'ColorScale', 'log'); caxis(h.axR, opts.clims_log);
    else
        caxis(h.axR, opts.clims);
    end
    cb = colorbar(h.axR); cb.Label.Interpreter = 'latex'; cb.Label.String = 'counts';
    title(h.axR, 'raw LIF emission', 'Interpreter', 'latex');
    [h.sfR, h.z9R, h.pnR] = overlays(h.axR, y_cm, nx, c99, cpen);

    % 2. concentration
    h.axC = nexttile(tl, 7, [1 5]);
    h.imC = imagesc(h.axC, y_cm, z_cm, nan(ny, nx));
    colormap(h.axC, cmapM);
    if opts.C_log
        set(h.axC, 'ColorScale', 'log'); caxis(h.axC, opts.C_clims_log);
    else
        caxis(h.axC, opts.C_clims);
    end
    cb = colorbar(h.axC); cb.Label.Interpreter = 'latex'; cb.Label.String = opts.C_label;
    title(h.axC, '$C$ (rhodamine, sheet profile removed)', 'Interpreter', 'latex');
    [h.sfC, h.z9C, h.pnC] = overlays(h.axC, y_cm, nx, c99, cpen);
    legend(h.axC, {'surface', '$z_{99}(y)$', 'dye penetration'}, ...
           'Interpreter', 'latex', 'Location', 'southeast', 'TextColor', 'k');

    for ax = [h.axR h.axC]
        set(ax, 'YDir', 'normal', 'Color', 'k');
        daspect(ax, [1 1 1]);                        % cm for cm, like axis image
        if isempty(opts.ylim_cm), xlim(ax, sort(y_cm([1 end]))); else, xlim(ax, opts.ylim_cm); end
        if isempty(opts.zlim_cm), ylim(ax, sort(z_cm([1 end]))); else, ylim(ax, opts.zlim_cm); end
        xlabel(ax, '$y$ (cm, cross-tank)', 'Interpreter', 'latex');
        ylabel(ax, '$z$ (cm, 0 = frame-1 water level)', 'Interpreter', 'latex');
        set(ax, 'FontSize', 12, 'FontName', 'times');
    end

    % 3. profile
    h.axP = nexttile(tl, 6, [2 1]);
    hold(h.axP, 'on');
    h.bandP = patch(h.axP, [0 1 1 0], [NaN NaN NaN NaN], c99, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    xline(h.axP, 0, ':');
    h.prP = plot(h.axP, NaN, NaN, 'k-', 'LineWidth', 1.2);
    h.z9P = plot(h.axP, [0 1], [NaN NaN], '-', 'Color', c99, 'LineWidth', 2);
    hold(h.axP, 'off'); grid(h.axP, 'on'); box(h.axP, 'on');
    if ~isempty(opts.profile_xlim), xlim(h.axP, opts.profile_xlim); end
    ylim(h.axP, opts.profile_zlim);
    xlabel(h.axP, ['$\overline{C}$ (all $y$)'], 'Interpreter', 'latex');
    ylabel(h.axP, '$z - \eta$ (cm, 0 = local surface)', 'Interpreter', 'latex');
    legend(h.axP, [h.z9P h.bandP], {'$z_{99}$', '$z_{99}(y)$ 5--95\%'}, ...
           'Interpreter', 'latex', 'Location', 'southeast');
    set(h.axP, 'FontSize', 12, 'FontName', 'times');
    h.ttP = title(h.axP, '', 'Interpreter', 'latex');

    % 4. hot-wire time series, vertical line at the frame time
    h.hwLine = [];
    if has_hw
        h.axH = nexttile(tl, 13, [1 6]);
        plot(h.axH, opts.hw.t, opts.hw.Y, 'LineWidth', 0.8);
        hold(h.axH, 'on');
        h.hwLine = xline(h.axH, opts.hw.t(1), 'r-', 'LineWidth', 2);
        hold(h.axH, 'off'); grid(h.axH, 'on'); box(h.axH, 'on');
        if isfield(opts.hw, 'xlim') && ~isempty(opts.hw.xlim), xlim(h.axH, opts.hw.xlim);
        else, xlim(h.axH, opts.hw.t([1 end])); end
        xlabel(h.axH, 'time since wind start (s)');
        ylabel(h.axH, 'velocity (m/s)');
        legend(h.axH, [opts.hw.names, {'frame'}], 'Interpreter', 'latex', 'Location', 'northwest');
        set(h.axH, 'FontSize', 12, 'FontName', 'times');
    end

    h.opts = opts;
    h.tt   = title(tl, '', 'Interpreter', 'latex', 'FontSize', 15);
end

function [hs, h9, hp] = overlays(ax, y_cm, nx, c99, cpen)
    hold(ax, 'on');
    hs = plot(ax, y_cm, nan(1, nx), 'w--', 'LineWidth', 1);
    h9 = plot(ax, NaN, NaN, '-',  'Color', c99,  'LineWidth', 2);
    hp = plot(ax, y_cm, nan(1, nx), ':', 'Color', cpen, 'LineWidth', 1.2);
    hold(ax, 'off');
end
