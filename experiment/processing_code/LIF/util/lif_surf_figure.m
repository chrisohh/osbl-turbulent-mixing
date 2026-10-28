function h = lif_surf_figure(x_cm, z_cm, opts)
%LIF_SURF_FIGURE Three-panel Water_SURF figure: raw | theta | profile.
%
%   h = LIF_SURF_FIGURE(x_cm, z_cm, opts)   builds the (empty) figure once
%   LIF_SURF_UPDATE(h, ...)                 fills it for one frame
%
% Shared by surf_frame_calibrated.m (one frame) and surf_theta_video.m
% (every frame), so the still and the video look the same.
%   1. transmission I/I_bg (opts.left = 'transmission', default: dye
%      dark, lighting removed, 1 = clear) or raw counts ('raw'), with
%      surface / z99(x)
%   2. theta = ln(I_bg/I), lighting removed (faint edge near z99 visible),
%      same overlays
%   3. horizontally averaged theta vs depth below the surface, all-x z99
%      and the 5-95% spread of the local z99(x)
%
% opts fields (defaults):
%   cmap 'bone', clims [350 600], theta_clims [0 0.05], theta_log false,
%   theta_clims_log [0.005 0.3], profile_xlim [-0.01 0.3],
%   profile_zlim [-15 0], fig_pos [40 40 2000 650], visible 'on',
%   left 'transmission' ('raw'), trans_clims [0.75 1.05],
%   hw [] -- hot-wire panel under the images (full width) with a vertical
%      line at the frame time (lif_surf_update's t_frame).  Struct with
%      t (s since wind start), Y ([nt x nseries] m/s), names (cell, LaTeX),
%      xlim ([] = whole record).  fig_pos grows to fit when given.
%   panels 'all' -- which panels to build:
%      'all'   the three panels above
%      'theta' theta only, full width, over the hot-wire panel.  Panels 1
%              and 3 are not created and their handles come back empty;
%              lif_surf_update skips whatever is missing.
%   show_z99 true -- false drops the z99 overlays (and the z99 entries from
%      the legends), for frames where z99 is not the point
%   wave [] -- wave-gauge eta panel under the hot-wire panel (same struct as
%      hw; lif_wind_clock's third output), same time axis and frame line
%   z99_band [5 95] -- percentiles of z99(x) shaded in the profile panel
%   ct [] -- stratified runs: lif_ct_profile output.  The pre-run CT profile
%      is drawn in the profile panel on its own top x axis, the pycnocline
%      shaded there and dashed on the theta image

    def = struct('cmap', 'bone', 'clims', [350 600], 'theta_clims', [0 0.05], ...
                 'theta_log', false, 'theta_clims_log', [0.005 0.3], ...
                 'profile_xlim', [-0.01 0.3], 'profile_zlim', [-15 0], ...
                 'fig_pos', [40 40 2000 650], 'visible', 'on', ...
                 'left', 'transmission', 'trans_clims', [0.75 1.05], 'hw', [], ...
                 'img_xlim', [], 'img_ylim', [], 'z99_floor_cm', [], ...
                 'panels', 'all', 'show_z99', true, 'ct', [], 'wave', [], 'z99_band', [5 95]);
    if nargin < 3, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end
    c99 = [0.3 0.75 1];
    nx = numel(x_cm);  ny = numel(z_cm);

    h.opts = opts;
    has_hw     = ~isempty(opts.hw);
    series     = {opts.hw, opts.wave};  series = series(~cellfun(@isempty, series));   % wind above wave gauge
    nser       = numel(series);
    theta_only = strcmp(opts.panels, 'theta');
    if nser > 0 && isequal(opts.fig_pos, def.fig_pos), opts.fig_pos(4) = 1000 + 250*(nser - 1); end
    h.fig  = figure('Position', opts.fig_pos, 'Color', 'w', 'Visible', opts.visible);
    ncol   = 7;
    img_rows = 1 + (nser > 0);                         % images: 2 rows over the series rows
    tl     = tiledlayout(h.fig, img_rows + nser, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');

    [h.axR, h.imR, h.sfR, h.z9R, h.bdR, h.unR, h.bdT, h.unT] = deal(gobjects(0));
    [h.axP, h.prP, h.z9P, h.bandP, h.ttP] = deal(gobjects(0));

    % 1. transmission I/I_bg (or raw counts)
    if ~theta_only
        h.axR = nexttile(tl, 1, [img_rows 3]);
        h.imR = imagesc(h.axR, x_cm, z_cm, nan(ny, nx));
        colormap(h.axR, opts.cmap);
        cb = colorbar(h.axR); cb.Label.Interpreter = 'latex';
        if strcmp(opts.left, 'transmission')
            caxis(h.axR, opts.trans_clims);
            cb.Label.String = '$I/I_{\rm bg}$';
            title(h.axR, 'transmission $I/I_{\rm bg}$ (lighting removed)', 'Interpreter', 'latex');
        else
            caxis(h.axR, opts.clims);
            cb.Label.String = 'counts';
            title(h.axR, 'raw', 'Interpreter', 'latex');
        end
        [h.sfR, h.z9R, h.bdR, h.unR] = overlays(h.axR, x_cm, nx, c99, opts.show_z99, ...
                                                strcmp(opts.left, 'transmission'));
    end

    % 2. theta -- full width when it is the only image
    if theta_only
        h.axT = nexttile(tl, 1, [img_rows ncol]);
    else
        h.axT = nexttile(tl, 4, [img_rows 3]);
    end
    h.imT = imagesc(h.axT, x_cm, z_cm, nan(ny, nx));
    colormap(h.axT, flipud(feval(opts.cmap, 256)));   % dye dark, as in the raw image
    if opts.theta_log
        set(h.axT, 'ColorScale', 'log'); caxis(h.axT, opts.theta_clims_log);
    else
        caxis(h.axT, opts.theta_clims);
    end
    cb = colorbar(h.axT); cb.Label.Interpreter = 'latex';  h.cbT = cb;
    cb.Label.String = '$\theta = \ln(I_{\rm bg}/I)$';
    [h.sfT, h.z9T, h.bdT, h.unT] = overlays(h.axT, x_cm, nx, c99, opts.show_z99, true);
    title(h.axT, '$\theta$ (dye concentration)', 'Interpreter', 'latex');
    if opts.show_z99
        legend(h.axT, [h.sfT h.z9T], {'surface', '$z_{99}(x)$'}, ...
               'Interpreter', 'latex', 'Location', 'southeast');
    else
        legend(h.axT, h.sfT, {'surface'}, 'Interpreter', 'latex', 'Location', 'southeast');
    end

    for ax = [h.axR h.axT]
        axis(ax, 'image'); set(ax, 'YDir', 'normal');
        if ax == h.axT || strcmp(opts.left, 'transmission')
            set(ax, 'Color', [0.78 0.78 0.78]);   % air above the surface (NaN) shows grey
        end
        % Requested limits, clipped to this run's image: runs with a smaller
        % field of view (run 92: 0.079 vs 0.125 mm/px) left empty margins
        % with limits set for the wide view
        xr = [min(x_cm) max(x_cm)];  zr = [min(z_cm) max(z_cm)];
        if ~isempty(opts.img_xlim)
            xlim(ax, [max(opts.img_xlim(1), xr(1)) min(opts.img_xlim(2), xr(2))]);
        else
            xlim(ax, xr);
        end
        if ~isempty(opts.img_ylim)
            ylim(ax, [max(opts.img_ylim(1), zr(1)) min(opts.img_ylim(2), zr(2))]);
        else
            ylim(ax, zr);
        end
        xlabel(ax, '$x$ (cm)', 'Interpreter', 'latex');
        ylabel(ax, '$z$ (cm, 0 = frame-1 water level)', 'Interpreter', 'latex');
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
        % from the wind start: the wave-gauge logger starts ~10 s before
        % the fan (hwT_wind from -10 s), the hot-wire series at 0
        xl = [max(0, min(cellfun(@(s) s.t(1), series))) max(cellfun(@(s) s.t(end), series))];
        if isfield(opts.hw, 'xlim') && ~isempty(opts.hw.xlim), xl = opts.hw.xlim; end
        linkaxes(axS, 'x');  xlim(axS(1), xl);
        xlabel(axS(end), 'time since wind start (s)');
        set(axS(1:end-1), 'XTickLabel', []);
        h.axH = axS(1);
    end

    % 3. profile: same height and z range as the theta image.  Its tile only
    % reserves the width; the axes is placed over it, at the theta image's
    % plot box (axis image shrinks that inside the tile, so a tile-filling
    % profile was taller than theta and on a different z scale)
    if ~theta_only
        hold_tile = nexttile(tl, 7, [img_rows 1]);  axis(hold_tile, 'off');
        % the figure title takes its height from the layout: it must exist
        % (with text) before measuring, or the images move down under the
        % already-placed profile once the first frame sets the title
        h.tt = title(tl, 'CoreView', 'Interpreter', 'latex', 'FontSize', 15);
        drawnow;
        fp = getpixelposition(h.fig);  fw = fp(3);  fh = fp(4);
        pt = getpixelposition(h.axT, true);
        xr = diff(xlim(h.axT));  zr = diff(ylim(h.axT));
        sc = min(pt(3) / xr, pt(4) / zr);              % px per cm in the theta image
        bh = zr * sc;  by = pt(2) + (pt(4) - bh) / 2;
        % width: from past the theta colour bar and its label (+ room for
        % this panel's own y label) to the figure edge -- the tile alone was
        % ~60 px wide and its text spilled out
        cbr = (h.cbT.Position(1) + h.cbT.Position(3)) * fw;
        x0  = cbr + 0.075*fw;  x1 = 0.985*fw;
        opts.profile_zlim = ylim(h.axT);
        h.axP = axes(h.fig, 'Units', 'pixels', 'Position', [x0 by max(x1 - x0, 60) bh]);
        hold(h.axP, 'on');
        h.bandP = patch(h.axP, opts.profile_xlim([1 2 2 1]), [NaN NaN NaN NaN], c99, ...
                        'FaceAlpha', 0.25, 'EdgeColor', 'none');
        xline(h.axP, 0, ':');
        h.prP = plot(h.axP, NaN, NaN, 'k-', 'LineWidth', 1.2);
        h.z9P = plot(h.axP, opts.profile_xlim, [NaN NaN], '-', 'Color', c99, 'LineWidth', 2);
        hold(h.axP, 'off'); grid(h.axP, 'on'); box(h.axP, 'on');
        xlim(h.axP, opts.profile_xlim); ylim(h.axP, opts.profile_zlim);
        xlabel(h.axP, '$\overline{\theta}$ (all $x$)', 'Interpreter', 'latex');
        ylabel(h.axP, '$z - \eta$ (cm, 0 = local surface)', 'Interpreter', 'latex');
        set(h.axP, 'FontSize', 12, 'FontName', 'times');
        % z99 value and legend go in the empty top of the panel (z - eta > 0
        % is air): a title above it spilled out of the narrow panel and
        % collided with the CT axis, and a legend lower down covered the profile
        h.ttP = text(h.axP, 0.04, 0.99, '', 'Units', 'normalized', 'Interpreter', 'latex', ...
                     'VerticalAlignment', 'top', 'FontSize', 11);
        lh = gobjects(0);  ls = {};
        if opts.show_z99
            lh = h.z9P;  ls = {'$z_{99}$'};
            if ~isempty(opts.z99_band)                 % [] = no spread band
                lh(end+1) = h.bandP;  ls{end+1} = sprintf('$z_{99}(x)$ %g--%g\\%%', opts.z99_band);
            end
        end

        % Stratified runs: pre-run CT profile on its own (top) x axis in the
        % same tile, pycnocline shaded on the profile and marked on theta
        if ~isempty(opts.ct)
            P = opts.ct;  cct = [0.85 0.33 0.10];
            hold(h.axP, 'on');
            patch(h.axP, opts.profile_xlim([1 2 2 1]), P.pyc([1 1 2 2]), cct, ...
                  'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            hold(h.axP, 'off');
            h.axCT = axes(h.fig, 'Units', 'pixels', 'Position', h.axP.Position);
            plot(h.axCT, P.v, P.z, '.-', 'Color', cct, 'LineWidth', 1, 'MarkerSize', 8);
            set(h.axCT, 'Color', 'none', 'XAxisLocation', 'top', 'YAxisLocation', 'right', ...
                'YTickLabel', [], 'XColor', cct, 'Box', 'off', 'FontSize', 10, 'FontName', 'times');
            ylim(h.axCT, opts.profile_zlim);
            linkprop([h.axP h.axCT], 'YLim');
            yline(h.axT, P.pyc, '--', 'Color', cct, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            hold(h.axP, 'on');                         % legend entry for the CT curve (top axis, V)
            lh(end+1) = plot(h.axP, NaN, NaN, '.-', 'Color', cct, 'MarkerSize', 8);
            ls{end+1} = 'CT (V, top)';
            hold(h.axP, 'off');
            uistack(h.axP, 'top');  set(h.axP, 'Color', 'none');   % legend and text above the CT axes
        end
        set(findall(h.fig, 'Type', 'axes', 'Units', 'pixels'), 'Units', 'normalized');   % scale with the figure
        if ~isempty(lh)
            lg = legend(h.axP, lh, ls, 'Interpreter', 'latex', 'FontSize', 8, 'Box', 'off', 'Location', 'northeast');
            drawnow;                                   % its size is only known once laid out
            lg.Location = 'none';
            ap = h.axP.Position;  lp = lg.Position;
            lg.Position = [ap(1) + ap(3) - lp(3), ap(2) + 0.885*ap(4) - lp(4), lp(3), lp(4)];   % under the z99 text
        end
    end

    h.opts = opts;

    if ~isfield(h, 'tt'), h.tt = title(tl, '', 'Interpreter', 'latex', 'FontSize', 15); end
end

function [hs, h9, hb, hu] = overlays(ax, x_cm, nx, c99, show_z99, grey)
    % grey layer (drawn first, under the lines) over the unlit columns.  The
    % air above the surface is grey through the axes colour (NaN pixels).
    hold(ax, 'on');
    [hb, hu] = deal(gobjects(0));
    if grey
        g  = [0.78 0.78 0.78];
        % (no surface-band layer: only the air is grey, by request; the
        % band is still left out of the z99 integral in lif_z99)
        hu = patch(ax, NaN, NaN, g, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hs = plot(ax, x_cm, nan(1, nx), 'c--', 'LineWidth', 1);
    h9 = gobjects(0);
    if show_z99
        h9 = plot(ax, NaN, NaN, '-', 'Color', c99, 'LineWidth', 2);
    end
    hold(ax, 'off');
end
