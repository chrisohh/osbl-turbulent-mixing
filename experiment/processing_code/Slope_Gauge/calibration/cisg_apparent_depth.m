function depth = cisg_apparent_depth(ref_image, cal, cfg)
%CISG_APPARENT_DEPTH Camera height and air-equivalent pattern depth from the colour target.
%
%   depth = CISG_APPARENT_DEPTH(ref_image, cal)
%   depth = CISG_APPARENT_DEPTH(ref_image, cal, cfg)
%
% Paraxial, camera looking down.  The surface-plane cal (cisg_plate_calibrate)
% gives mm/px at the water surface, so the camera height above the surface is
%
%   m = pitch / mm_per_px,   H = f (1 + 1/m).
%
% The colour target (cfg.target_mm square) sits deeper than the surface, so it
% looks smaller: its left and right edges, taken to the surface plane with the
% cal, are W_s apart.  With r = target_mm / W_s,
%
%   L_eff = H (r - 1) = h/n_w + t/n_g + g          (air-equivalent depth)
%   h_est = n_w (L_eff - t/n_g - g)                (water depth, sanity check)
%
% (t, n_g: glass bottom; g: air gap to the target).  The slopes use L_eff and the
% pattern-plane pixel size (mm/px at the surface * r), not t or g.
%
% Edges: the target is a smooth colour ramp on a dark surround, blurred over
% ~100 px, so a hard threshold fails.  For each row of the middle band the
% edge is the strongest relative (log gray) gradient in an outer window; a
% robust line col(row) is fitted per side (handles a slightly rotated target),
% both lines go to the surface plane and W_s is their separation along the
% mean normal.
%
% ref_image - RGB reference (still-water) image
% cal       - cisg_plate_calibrate output
% cfg fields (defaults):
%   f_mm 35, pixel_um 5.5, target_mm 508, n_w 1.333
%   t_mm 19, n_g 1.52, gap_mm 76.2     glass bottom and air gap to the target
%   H_mm []             camera height override (measured); [] = from the cal
%   edge_cols []        [left right] column override at the centre row
%   row_band [0.25 0.75]  fraction of image rows used
%   win_left [0.04 0.12], win_right [0.88 0.99]   edge search windows (fraction of width)
%   smooth_px [4 15]    gaussian sigma [rows cols] before the gradient
%   fig_path ''         save a diagnostic figure here (png), '' = none

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(fileparts(here)), 'IR Camera', 'calibration'));   % surf_px2world

    def = struct('f_mm', 35, 'pixel_um', 5.5, 'target_mm', 508, 'n_w', 1.333, ...
                 't_mm', 19, 'n_g', 1.52, 'gap_mm', 76.2, 'H_mm', [], 'edge_cols', [], ...
                 'row_band', [0.25 0.75], 'win_left', [0.04 0.12], 'win_right', [0.88 0.99], ...
                 'smooth_px', [4 15], 'fig_path', '');
    if nargin < 3, cfg = struct(); end
    for f = fieldnames(def).'
        if ~isfield(cfg, f{1}), cfg.(f{1}) = def.(f{1}); end
    end

    %% Camera height
    mmpp = cal.mm_per_px_center;                       % surface mm/px at the plate (~ image) centre
    m = (cfg.pixel_um / 1000) / mmpp;
    H_cal = cfg.f_mm * (1 + 1 / m);
    if isempty(cfg.H_mm), H = H_cal; H_src = 'cal'; else, H = cfg.H_mm; H_src = 'cfg'; end

    %% Target edges
    [nr, nc, ~] = size(ref_image);
    rows = round(cfg.row_band(1) * nr) : round(cfg.row_band(2) * nr);
    rows = rows(1:8:end);
    g = mean(double(ref_image(rows(1):rows(end), :, :)), 3);
    g = imgaussfilt(g, cfg.smooth_px);
    d = gradient(log(g + 5));                          % relative gradient along columns
    d = d(rows - rows(1) + 1, :);
    wl = round(cfg.win_left * nc);   wl = wl(1):wl(2);
    wr = round(cfg.win_right * nc);  wr = wr(1):wr(2);
    [~, il] = max(d(:, wl), [], 2);   cl = wl(il)';    % rising edge on the left
    [~, ir] = min(d(:, wr), [], 2);   cr = wr(ir)';    % falling edge on the right
    rc = (nr + 1) / 2;
    [pl, okl] = robust_line(rows(:), cl);
    [pr, okr] = robust_line(rows(:), cr);
    ec = [polyval(pl, rc) polyval(pr, rc)];            % edge columns at the centre row
    if ~isempty(cfg.edge_cols), ec = cfg.edge_cols; pl(2) = ec(1) - pl(1)*rc; pr(2) = ec(2) - pr(1)*rc; end

    % to the surface plane: points on each edge line at centre row +/- 500 rows
    dr = 500;
    line_uv = @(p) [polyval(p, [rc - dr; rc + dr]), [rc - dr; rc + dr]];
    xyL = surf_px2world(cal, line_uv(pl));
    xyR = surf_px2world(cal, line_uv(pr));
    eL = diff(xyL) / norm(diff(xyL));  eR = diff(xyR) / norm(diff(xyR));
    e = eL + eR;  e = e / norm(e);                     % mean edge direction
    nrm = [e(2) -e(1)];                                % normal (tank x ~ image right)
    pL = surf_px2world(cal, [ec(1) rc]);  pR = surf_px2world(cal, [ec(2) rc]);
    W_s = abs((pR - pL) * nrm');
    tilt = atand(e(1) / e(2));                         % edge angle from the tank y axis (deg)

    %% Depth
    r = cfg.target_mm / W_s;
    L_eff = H * (r - 1);
    h_est = cfg.n_w * (L_eff - cfg.t_mm / cfg.n_g - cfg.gap_mm);

    depth = struct('H_mm', H, 'H_source', H_src, 'H_from_cal_mm', H_cal, 'm', m, ...
        'mm_per_px_surface', mmpp, 'W_s_mm', W_s, 'edge_cols', ec, ...
        'edge_cols_rows', [cl cr], 'edge_rows', rows(:), 'edge_line_px', [pl; pr], ...
        'edge_inlier_frac', [mean(okl) mean(okr)], 'edge_tilt_deg', tilt, ...
        'r', r, 'L_eff_mm', L_eff, 'h_est_mm', h_est, ...
        'pattern_mm_per_px', mmpp * r, 'cfg', cfg, ...
        'h_est_note', sprintf('h_est depends on ASSUMED glass t = %.1f mm (n_g %.2f) and air gap %.1f mm; slopes use L_eff only', ...
                              cfg.t_mm, cfg.n_g, cfg.gap_mm));

    fprintf('Apparent depth: H %.0f mm (%s; cal %.0f), m %.4f, edges %.0f / %.0f px (inliers %.0f%% / %.0f%%, tilt %.2f deg)\n', ...
            H, H_src, H_cal, m, ec, 100*depth.edge_inlier_frac, tilt);
    fprintf('  W_s %.1f mm -> r %.4f; L_eff %.1f mm (used for slopes)\n', W_s, r, L_eff);
    fprintf('  h_est %.1f mm: sanity check only, depends on ASSUMED t = %.1f mm, g = %.1f mm\n', ...
            h_est, cfg.t_mm, cfg.gap_mm);

    %% Diagnostic figure
    if ~isempty(cfg.fig_path)
        fig = figure('Visible', 'off', 'Position', [100 100 1200 900]);
        subplot(2,1,1);
        imagesc(uint8(ref_image)); axis image; hold on;
        plot(cl, rows, 'g.', polyval(pl, [1 nr]), [1 nr], 'g-', 'LineWidth', 1);
        plot(cr, rows, 'r.', polyval(pr, [1 nr]), [1 nr], 'r-', 'LineWidth', 1);
        title(sprintf('Reference image: target edges %.0f / %.0f px, W_s %.1f mm, r %.3f, L_{eff} %.0f mm', ec, W_s, r, L_eff));
        subplot(2,1,2);
        pr_ = squeeze(mean(double(ref_image(rows(1):rows(end), :, :)), 1));
        plot(pr_(:,1), 'r'); hold on; plot(pr_(:,2), 'g'); plot(pr_(:,3), 'b');
        xline(ec(1), 'k--'); xline(ec(2), 'k--');
        xlabel('column (px)'); ylabel('mean count over the row band'); grid on;
        exportgraphics(fig, cfg.fig_path, 'Resolution', 100);
        close(fig);
    end
end

function [p, ok] = robust_line(x, y)
    % col = p(1) * row + p(2), dropping outliers (3 MAD) twice
    ok = true(size(y));
    for it = 1:3
        p = polyfit(x(ok), y(ok), 1);
        r = y - polyval(p, x);
        s = max(1.4826 * median(abs(r(ok) - median(r(ok)))), 2);
        ok = abs(r) <= 3 * s;
    end
    p = polyfit(x(ok), y(ok), 1);
end
