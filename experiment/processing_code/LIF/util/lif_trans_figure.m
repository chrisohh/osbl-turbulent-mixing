function h = lif_trans_figure(y_cm, z_cm, opts)
%LIF_TRANS_FIGURE Three-panel transverse-LIF figure: raw | C | profile.
%
%   h = LIF_TRANS_FIGURE(y_cm, z_cm, opts)   builds the (empty) figure once
%   LIF_TRANS_UPDATE(h, ...)                 fills it for one frame
%
% The transverse (y, z) counterpart of lif_surf_figure, laid out the same
% way so the two videos read alike: raw | C side by side, the profile on
% the right at the C image's height and z range, wind and wave-gauge panels
% underneath, CT profile and pycnocline for stratified runs.  Shared by
% trans_frame_calibrated.m (one frame) and trans_conc_video.m (every frame).
% Differences from the side view, all because the dye emits rather than
% absorbs:
%   * the colour map is NOT flipped -- dye is bright, like the raw image,
%     and the background (air, no dye) is black;
%   * the raw panel is log-scaled by default, since the surface line is
%     ~150 counts and the deep plumes ~1 count over a black background;
%   * the overlays are the surface and z99(y), plus the dye-penetration
%     contour z_pen(y) from lif_trans_conc.
%
%   1. raw counts (opts.raw_log: log colour scale)
%   2. concentration C from lif_trans_conc, in the units set by its
%      normalize option
%   3. y-averaged C vs depth below the local surface, z99 and the spread of
%      z99(y)
%
% opts fields (defaults):
%   cmap 'inferno' (falls back to 'hot' if not on the path), clims [0 150],
%   raw_log true, clims_log [0.5 300], C_clims [0 1], C_log false,
%   C_clims_log [1e-3 1], C_label '$C$ (normalised)',
%   profile_xlim [] ([] = autoscale per frame in lif_trans_update),
%   fig_pos [40 40 2000 650] (grows with the time-series panels),
%   visible 'on', show_ds 1 (display every show_ds-th pixel),
%   zlim_cm [] / ylim_cm [] -- image-panel limits ([] = whole frame).  The
%      analysis wants the WHOLE frame (lif_trans_conc needs clear water
%      well below the dye for its noise floor); zoom here, not in the data.
%      The profile takes the image's z range.
%   hw [], wave [] -- wind and wave-gauge panels (lif_wind_clock outputs)
%      with a vertical line at the frame time, as in lif_surf_figure
%   ct [] -- stratified runs: lif_ct_profile output (top x axis of the
%      profile panel, pycnocline shaded there and dashed on the C image)
%   z99_band [5 95] -- percentiles of z99(y) shaded in the profile ([] = none)
%   z99_floor_cm [] -- per-frame "resolved" gate (lif_z99_gate), set by
%      the caller in h.opts before each lif_trans_update

    def = struct('cmap', 'inferno', 'clims', [0 150], ...
                 'raw_log', true, 'clims_log', [0.5 300], ...
                 'C_clims', [0 1], 'C_log', false, 'C_clims_log', [1e-3 1], ...
                 'C_label', '$C$ (normalised)', ...
                 'profile_xlim', [], 'profile_zlim', [], ...
                 'fig_pos', [40 40 2000 650], 'visible', 'on', ...
                 'zlim_cm', [], 'ylim_cm', [], 'hw', [], 'wave', [], 'ct', [], ...
                 'z99_band', [5 95], 'z99_floor_cm', [], 'show_ds', 1);
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

    series = {opts.hw, opts.wave};  series = series(~cellfun(@isempty, series));   % wind above wave gauge
    nser   = numel(series);
    % the transverse images are ~2:1 (wider than the side view's 1.6:1), so
    % less height than lif_surf_figure, or they float in white space
    if nser > 0 && isequal(opts.fig_pos, def.fig_pos), opts.fig_pos(4) = 500 + 250*nser; end
    % painters: with OpenGL a video frame took 1 s or ~50 s to print
    % depending on which OpenGL driver MATLAB picked at start-up on this
    % machine; painters avoids OpenGL altogether (~1.5 s, every session)
    h.fig = figure('Position', opts.fig_pos, 'Color', 'w', 'Visible', opts.visible, ...
                   'Renderer', 'painters');
    ncol     = 7;
    img_rows = 1 + (nser > 0);                         % images: 2 rows over the series rows
    tl = tiledlayout(h.fig, img_rows + nser, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 1. raw counts
    h.axR = nexttile(tl, 1, [img_rows 3]);
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
    h.axC = nexttile(tl, 4, [img_rows 3]);
    h.imC = imagesc(h.axC, y_cm, z_cm, nan(ny, nx));
    colormap(h.axC, cmapM);
    if opts.C_log
        set(h.axC, 'ColorScale', 'log'); caxis(h.axC, opts.C_clims_log);
    else
        caxis(h.axC, opts.C_clims);
    end
    cb = colorbar(h.axC); cb.Label.Interpreter = 'latex'; cb.Label.String = opts.C_label;  h.cbC = cb;
    title(h.axC, '$C$ (dye concentration, sheet profile removed)', 'Interpreter', 'latex');
    [h.sfC, h.z9C, h.pnC] = overlays(h.axC, y_cm, nx, c99, cpen);
    legend(h.axC, [h.sfC h.z9C h.pnC], {'surface', '$z_{99}(y)$', 'dye penetration'}, ...
           'Interpreter', 'latex', 'Location', 'southeast', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'w');

    for ax = [h.axR h.axC]
        set(ax, 'YDir', 'normal', 'Color', 'k');
        daspect(ax, [1 1 1]);                        % cm for cm, like axis image
        if isempty(opts.ylim_cm), xlim(ax, sort(y_cm([1 end]))); else, xlim(ax, opts.ylim_cm); end
        if isempty(opts.zlim_cm), ylim(ax, sort(z_cm([1 end]))); else, ylim(ax, opts.zlim_cm); end
        xlabel(ax, '$y$ (cm, cross-tank)', 'Interpreter', 'latex');
        ylabel(ax, '$z$ (cm, 0 = reference water level)', 'Interpreter', 'latex');
        set(ax, 'FontSize', 12, 'FontName', 'times');
    end

    % 4. wind and wave-gauge time series, vertical line at the frame time
    %    (built before the profile: the profile is sized to the final layout)
    h.hwLine = gobjects(0);  axS = gobjects(0);
    for k = 1:nser
        s = series{k};
        ax = nexttile(tl, ncol*(img_rows + k - 1) + 1, [1 ncol]);
        plot(ax, s.t, s.Y, 'LineWidth', 0.8);
        hold(ax, 'on');
        h.hwLine(end+1) = xline(ax, s.t(1), 'r-', 'LineWidth', 2, 'HandleVisibility', 'off');   % moved per frame
        hold(ax, 'off'); grid(ax, 'on'); box(ax, 'on');
        if isfield(s, 'ylabel'), ylabel(ax, s.ylabel); else, ylabel(ax, 'velocity (m/s)'); end
        legend(ax, s.names, 'Interpreter', 'latex', 'Location', 'northwest');
        set(ax, 'FontSize', 12, 'FontName', 'times');
        axS(end+1) = ax; %#ok<AGROW>
    end
    if nser > 0
        % from the wind start: the wave-gauge logger starts ~10 s before the fan
        xl = [max(0, min(cellfun(@(s) s.t(1), series))) max(cellfun(@(s) s.t(end), series))];
        if isfield(opts.hw, 'xlim') && ~isempty(opts.hw.xlim), xl = opts.hw.xlim; end
        linkaxes(axS, 'x');  xlim(axS(1), xl);
        xlabel(axS(end), 'time since wind start (s)');
        set(axS(1:end-1), 'XTickLabel', []);
        h.axH = axS(1);
    end

    % 3. profile: same height and z range as the C image.  Its tile only
    % reserves the width; the axes is placed over it at the C image's plot
    % box (daspect shrinks that inside the tile).  The figure title must
    % exist (with text) before measuring, or the images move down under the
    % already-placed profile once the first frame sets it.
    hold_tile = nexttile(tl, 7, [img_rows 1]);  axis(hold_tile, 'off');
    h.tt = title(tl, 'CoreView', 'Interpreter', 'latex', 'FontSize', 15);
    drawnow;
    fp = getpixelposition(h.fig);  fw = fp(3);
    pt = getpixelposition(h.axC, true);
    yr = diff(xlim(h.axC));  zr = diff(ylim(h.axC));
    sc = min(pt(3) / yr, pt(4) / zr);                  % px per cm in the C image
    bh = zr * sc;  by = pt(2) + (pt(4) - bh) / 2;
    cbr = (h.cbC.Position(1) + h.cbC.Position(3)) * fw;
    x0  = cbr + 0.075*fw;  x1 = 0.985*fw;
    opts.profile_zlim = ylim(h.axC);
    h.axP = axes(h.fig, 'Units', 'pixels', 'Position', [x0 by max(x1 - x0, 60) bh]);
    px = opts.profile_xlim;  if isempty(px), px = [0 1]; end
    hold(h.axP, 'on');
    h.bandP = patch(h.axP, px([1 2 2 1]), [NaN NaN NaN NaN], c99, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    xline(h.axP, 0, ':');
    h.prP = plot(h.axP, NaN, NaN, 'k-', 'LineWidth', 1.2);
    h.z9P = plot(h.axP, px, [NaN NaN], '-', 'Color', c99, 'LineWidth', 2);
    hold(h.axP, 'off'); grid(h.axP, 'on'); box(h.axP, 'on');
    xlim(h.axP, px);  ylim(h.axP, opts.profile_zlim);
    xlabel(h.axP, '$\overline{C}$ (all $y$)', 'Interpreter', 'latex');
    ylabel(h.axP, '$z - \eta$ (cm, 0 = local surface)', 'Interpreter', 'latex');
    set(h.axP, 'FontSize', 12, 'FontName', 'times');
    % z99 value and legend in the empty top of the panel (z - eta > 0 is air)
    h.ttP = text(h.axP, 0.04, 0.99, '', 'Units', 'normalized', 'Interpreter', 'latex', ...
                 'VerticalAlignment', 'top', 'FontSize', 11);
    lh = h.z9P;  ls = {'$z_{99}$'};
    if ~isempty(opts.z99_band)
        lh(end+1) = h.bandP;  ls{end+1} = sprintf('$z_{99}(y)$ %g--%g\\%%', opts.z99_band);
    end

    % Stratified runs: pre-run CT profile on its own (top) x axis in the
    % same box, pycnocline shaded on the profile and dashed on the C image
    h.axCT = gobjects(0);
    if ~isempty(opts.ct)
        P = opts.ct;  cct = [0.85 0.33 0.10];
        hold(h.axP, 'on');
        h.pycP = patch(h.axP, px([1 2 2 1]), P.pyc([1 1 2 2]), cct, ...
                       'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        hold(h.axP, 'off');
        h.axCT = axes(h.fig, 'Units', 'pixels', 'Position', h.axP.Position);
        plot(h.axCT, P.v, P.z, '.-', 'Color', cct, 'LineWidth', 1, 'MarkerSize', 8);
        set(h.axCT, 'Color', 'none', 'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
            'YTickLabel', [], 'XColor', cct, 'Box', 'off', 'FontSize', 10, 'FontName', 'times');
        ylim(h.axCT, opts.profile_zlim);
        linkprop([h.axP h.axCT], 'YLim');
        yline(h.axC, P.pyc, '--', 'Color', cct, 'LineWidth', 0.8, 'HandleVisibility', 'off');
        hold(h.axP, 'on');                             % legend entry for the CT curve (top axis, V)
        lh(end+1) = plot(h.axP, NaN, NaN, '.-', 'Color', cct, 'MarkerSize', 8);
        ls{end+1} = 'CT (V, top)';
        hold(h.axP, 'off');
        uistack(h.axP, 'top');  set(h.axP, 'Color', 'none');   % legend and text above the CT axes
    end
    set(findall(h.fig, 'Type', 'axes', 'Units', 'pixels'), 'Units', 'normalized');   % scale with the figure
    lg = legend(h.axP, lh, ls, 'Interpreter', 'latex', 'FontSize', 8, 'Box', 'off', 'Location', 'northeast');
    drawnow;                                           % its size is only known once laid out
    lg.Location = 'none';
    ap = h.axP.Position;  lp = lg.Position;
    lg.Position = [ap(1) + ap(3) - lp(3), ap(2) + 0.885*ap(4) - lp(4), lp(3), lp(4)];   % under the z99 text

    h.opts = opts;
end

function [hs, h9, hp] = overlays(ax, y_cm, nx, c99, cpen)
    hold(ax, 'on');
    hs = plot(ax, y_cm, nan(1, nx), 'c--', 'LineWidth', 1);
    h9 = plot(ax, NaN, NaN, '-',  'Color', c99,  'LineWidth', 2);
    hp = plot(ax, y_cm, nan(1, nx), ':', 'Color', cpen, 'LineWidth', 1.2);
    hold(ax, 'off');
end
