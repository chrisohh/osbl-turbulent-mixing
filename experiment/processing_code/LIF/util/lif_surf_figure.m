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

    def = struct('cmap', 'bone', 'clims', [350 600], 'theta_clims', [0 0.05], ...
                 'theta_log', false, 'theta_clims_log', [0.005 0.3], ...
                 'profile_xlim', [-0.01 0.3], 'profile_zlim', [-15 0], ...
                 'fig_pos', [40 40 2000 650], 'visible', 'on', ...
                 'left', 'transmission', 'trans_clims', [0.75 1.05], 'hw', []);
    if nargin < 3, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end
    c99 = [0.3 0.75 1];
    nx = numel(x_cm);  ny = numel(z_cm);

    h.opts = opts;
    has_hw = ~isempty(opts.hw);
    if has_hw && isequal(opts.fig_pos, def.fig_pos), opts.fig_pos(4) = 1000; end
    h.fig  = figure('Position', opts.fig_pos, 'Color', 'w', 'Visible', opts.visible);
    tl     = tiledlayout(h.fig, 1 + 2*has_hw, 7, 'TileSpacing', 'compact', 'Padding', 'compact');
    img_rows = 1 + has_hw;                             % images get 2/3 of the height with hw

    % 1. transmission I/I_bg (or raw counts)
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
    [h.sfR, h.z9R] = overlays(h.axR, x_cm, nx, c99);

    % 2. theta
    h.axT = nexttile(tl, 4, [img_rows 3]);
    h.imT = imagesc(h.axT, x_cm, z_cm, nan(ny, nx));
    colormap(h.axT, flipud(feval(opts.cmap, 256)));   % dye dark, as in the raw image
    if opts.theta_log
        set(h.axT, 'ColorScale', 'log'); caxis(h.axT, opts.theta_clims_log);
    else
        caxis(h.axT, opts.theta_clims);
    end
    cb = colorbar(h.axT); cb.Label.Interpreter = 'latex';
    cb.Label.String = '$\theta = \ln(I_{\rm bg}/I)$';
    [h.sfT, h.z9T] = overlays(h.axT, x_cm, nx, c99);
    title(h.axT, '$\theta$ (dye concentration)', 'Interpreter', 'latex');
    legend(h.axT, {'surface', '$z_{99}(x)$'}, ...
           'Interpreter', 'latex', 'Location', 'southeast');

    for ax = [h.axR h.axT]
        axis(ax, 'image'); set(ax, 'YDir', 'normal');
        xlabel(ax, '$x$ (cm)', 'Interpreter', 'latex');
        ylabel(ax, '$z$ (cm, 0 = frame-1 water level)', 'Interpreter', 'latex');
        set(ax, 'FontSize', 12, 'FontName', 'times');
    end

    % 3. profile
    h.axP = nexttile(tl, 7, [img_rows 1]);
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
    legend(h.axP, [h.z9P h.bandP], {'$z_{99}$', '$z_{99}(x)$ 5--95\%'}, ...
           'Interpreter', 'latex', 'Location', 'southeast');
    set(h.axP, 'FontSize', 12, 'FontName', 'times');
    h.ttP = title(h.axP, '', 'Interpreter', 'latex');

    % 4. hot-wire time series, vertical line at the frame time
    h.hwLine = [];
    if has_hw
        h.axH = nexttile(tl, 7*img_rows + 1, [1 7]);
        plot(h.axH, opts.hw.t, opts.hw.Y, 'LineWidth', 0.8);
        hold(h.axH, 'on');
        h.hwLine = xline(h.axH, opts.hw.t(1), 'r-', 'LineWidth', 2);   % moved per frame
        hold(h.axH, 'off'); grid(h.axH, 'on'); box(h.axH, 'on');
        if isfield(opts.hw, 'xlim') && ~isempty(opts.hw.xlim), xlim(h.axH, opts.hw.xlim);
        else, xlim(h.axH, opts.hw.t([1 end])); end
        xlabel(h.axH, 'time since wind start (s)');
        ylabel(h.axH, 'velocity (m/s)');
        legend(h.axH, [opts.hw.names, {'frame'}], 'Interpreter', 'latex', 'Location', 'northwest');
        set(h.axH, 'FontSize', 12, 'FontName', 'times');
    end
    h.opts = opts;

    h.tt = title(tl, '', 'Interpreter', 'latex', 'FontSize', 15);
end

function [hs, h9] = overlays(ax, x_cm, nx, c99)
    hold(ax, 'on');
    hs = plot(ax, x_cm, nan(1, nx), 'c--', 'LineWidth', 1);
    h9 = plot(ax, NaN, NaN, '-', 'Color', c99, 'LineWidth', 2);
    hold(ax, 'off');
end
