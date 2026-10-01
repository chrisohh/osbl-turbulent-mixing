function cal = calibrate_surface_single(snap_path, opts)
%CALIBRATE_SURFACE_SINGLE  Pixel <-> water-surface map from ONE view of the
% full-size coarse board (calib_target_params('coarse'), 40 mm squares)
% floating on the water surface, e.g. Snap-000096.ats.
%
%   cal = calibrate_surface_single('D:\HLAB_2026\IR_camera\Snap-000096.ats');
%   cal = calibrate_surface_single(path, struct('fit_center', true));
%
% A single planar view cannot support the full Zhang calibration in
% calibrate_intrinsic.m, but it does not need to: everything we image lies
% on (approximately) the same plane as the board. The model is a planar
% homography plus an n-term radial lens distortion about a fixed centre,
%
%   p_u = H * [x; y; 1]                                (mm -> undistorted px)
%   p_d = c + (p_u - c) .* (1 + k1 r^2 + k2 r^4 + k3 r^6),  r = |p_u - c| / f
%
% fit by minimising corner reprojection error in raw pixels. H is solved
% linearly (fitgeotrans) inside the objective, so the nonlinear search is
% only over [k1..kn] (or [k1..kn cx cy] with opts.fit_center). The radial
% term is what replaces the side-wall polyfits of distortion_correction.m:
% it straightens the walls, and the homography handles the perspective
% (trapezoid) scaling that the wall method could only approximate with a
% uniform pixels-per-metre.
%
% Tank coordinates (mm): x = along-wind, y = cross-wind, on the surface.
%   y = 0 at the board centre = tank centreline (the 482.6 mm plate
%         self-centres in the 500 mm tank, so |error| < ~9 mm),
%         +y toward image right (same sign as distortion_correction.m).
%   x = 0 at the board centre (shift with opts.along_offset_mm),
%         +x toward image bottom (flip with opts.along_sign = -1).
% Board-to-tank axis assignment is chosen automatically from the image
% orientation, so it does not matter which way the square plate was laid.
%
% Options (all optional):
%   frame_idx        1        frame within the snapshot file
%   focal_px         520      13 mm / 25 um pitch; only normalises r
%   fit_center       false    also fit the distortion centre
%   n_radial         2        radial coefficients k1..kn. On Snap-000096,
%                             3 terms gained 0.006 px on the corners but set
%                             k3 = -1.2, which hooks the map inward and folds
%                             it beyond the last corner row; 1 term folds
%                             before the image corners (strong barrel).
%   along_offset_mm  0        tank x of the board centre
%   cross_offset_mm  0        tank y of the board centre
%   along_sign       +1       +1: x grows down the image, -1: up
%   wall_y_mm        250      side walls, for the overlay check
%   flat             []       board-free frame (same size) subtracted before
%                             detection, to remove the narcissus disc
%   use_lines        true     add hand-picked wall / pattern-edge lines
%   picks_path       calibration/line_picks_<snap name>.mat (reused if present)
%   repick           false    true = click all lines again; or a cell list
%                             of line names to redo only those, e.g.
%                             {'pattern bottom edge', 'left wall'}
%   wall_image       []       image to click the walls on ([] = snapshot);
%                             a recording frame shows the waterline better
%   line_weight      1        weight of a picked point vs a corner (px/px)
%   auto_walls       true     pre-detect the waterline (strongest cross-wise
%                             gradient near the predicted wall); Enter accepts
%   wall_band_px     25       search half-width around the predicted wall
%   out_path         calibration/surface_cal_<snap name>.mat ('' = no save)
%   show             true
%   square_mm        []       MEASURED corner-to-corner pitch; [] = design
%                             (calib_target_params). Measure across many
%                             squares (e.g. 10 = 400 mm design) or caliper the
%                             100 mm scale bar - one square's edges are biased
%                             by ablation spread, the corner pitch is not.

here = fileparts(mfilename('fullpath'));
addpath(here);

if nargin < 2, opts = struct(); end
[~, snap_name] = fileparts(snap_path);
def = struct('frame_idx', 1, 'focal_px', 520, 'fit_center', false, 'n_radial', 2, ...
    'along_offset_mm', 0, 'cross_offset_mm', 0, 'along_sign', 1, ...
    'wall_y_mm', 250, 'flat', [], ...
    'use_lines', true, 'picks_path', fullfile(here, ['line_picks_' snap_name '.mat']), ...
    'repick', false, 'wall_image', [], 'line_weight', 1, ...
    'auto_walls', true, 'wall_band_px', 25, ...
    'out_path', fullfile(here, ['surface_cal_' snap_name '.mat']), ...
    'show', true, 'square_mm', []);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end
end
opts.snap_path = snap_path;

target = calib_target_params('coarse');

%% ---- Load + detect ----
frame = read_calib_frame(snap_path, opts.frame_idx);
image_size = size(frame);
[pts, bs] = detect_board(frame, [target.n_rows target.n_cols], opts.flat);
fprintf('Detected %d corners (board %dx%d squares) in %s\n', ...
    size(pts,1), bs(1), bs(2), snap_name);

% Board-local corner coordinates, centred on the board
if isempty(opts.square_mm), opts.square_mm = target.square_mm; end
fprintf('Corner pitch %.3f mm (design %.1f mm)\n', opts.square_mm, target.square_mm);
wp = generateCheckerboardPoints(bs, opts.square_mm);
wp = wp - mean(wp, 1);

%% ---- Fit [k1..kn (cx cy)] with H solved linearly inside ----
cal0.center     = ([image_size(2) image_size(1)] + 1) / 2;
cal0.focal_px   = opts.focal_px;
cal0.k          = zeros(1, opts.n_radial);
cal0.T          = eye(3);
cal0.image_size = image_size;

rms_h = reproj_rms(cal0.k, pts, wp, cal0, false);

theta0 = [-0.05 zeros(1, opts.n_radial - 1)];
if opts.fit_center, theta0 = [theta0 cal0.center]; end
theta = fminsearch(@(th) reproj_rms(th, pts, wp, cal0, opts.fit_center), ...
    theta0, optimset('TolX', 1e-9, 'TolFun', 1e-9, ...
    'MaxFunEvals', 5000, 'MaxIter', 5000));

cal = cal0;
cal.k = theta(1:opts.n_radial);
if opts.fit_center, cal.center = theta(end-1:end); end

%% ---- Assign board axes to tank axes from the image orientation ----
cal = refit_T(cal, pts, wp);
p0 = surf_world2px(cal, [0 0]);
J  = [(surf_world2px(cal, [1 0]) - p0)', (surf_world2px(cal, [0 1]) - p0)'];

best = -Inf;
for P = {[1 0; 0 1], [0 1; 1 0]}
    for sgn = [1 1; 1 -1; -1 1; -1 -1]'
        S = diag(sgn) * P{1};                  % tank = S * board
        Jt = J * S';                            % columns: d(u,v)/dx, d(u,v)/dy
        score = opts.along_sign * Jt(2,1) + Jt(1,2);
        if score > best, best = score; S_best = S; end
    end
end
wp_tank = (S_best * wp')' + [opts.along_offset_mm opts.cross_offset_mm];
cal = refit_T(cal, pts, wp_tank);
rms_corners_only = reproj_rms(cal.k, pts, wp_tank, cal, false);

%% ---- Line constraints: side walls + outer pattern edges ----
% The 90 inner corners stop ~1 square short of the pattern edge and ~50 mm
% short of the walls, so on their own they leave the radial polynomial free
% to hook inward at the bottom corners. Hand-picked points on lines of
% known tank position pin that down: the waterline on each wall must map
% to y = -/+wall_y_mm (500 mm apart) and the outer edge of the checker
% pattern to its known x / y. H stays a linear fit to the corners; the
% lines enter through the radial terms (and centre, if fitted).
picks = [];
if opts.use_lines
    lines = define_lines(wp_tank, opts.square_mm, opts);
    if isempty(opts.wall_image), opts.wall_image = frame; end
    picks = get_line_picks(lines, frame, opts, cal);
end
if ~isempty(picks)
    th0 = cal.k;
    if opts.fit_center, th0 = [th0 cal.center]; end
    theta = fminsearch(@(th) joint_cost(th, pts, wp_tank, cal, opts.fit_center, ...
        picks, opts.line_weight), th0, optimset('TolX', 1e-9, 'TolFun', 1e-9, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000));
    cal.k = theta(1:opts.n_radial);
    if opts.fit_center, cal.center = theta(end-1:end); end
    cal = refit_T(cal, pts, wp_tank);
end
cal.line_picks = picks;

res = surf_world2px(cal, wp_tank) - pts;
cal.rms_px      = sqrt(mean(sum(res.^2, 2)));
cal.rms_px_homography_only = rms_h;
cal.board_to_tank = S_best;
cal.image_points  = pts;
cal.world_points  = wp_tank;
cal.board_size    = bs;
cal.snap_path     = snap_path;
cal.opts          = opts;

% World footprint of the full image (for choosing a rectification grid)
[uu, vv] = border_pixels(image_size);
xy_border = surf_px2world(cal, [uu vv]);
cal.border_mm = xy_border;
% Trusted footprint = detected corners + 60 mm (1.5 squares). Image-edge
% pixels are unconstrained: the top rows sit near the homography horizon
% and the radial polynomial extrapolates past the last corners, so the
% border maps to x ranges of hundreds of metres.
cal.bbox_mm   = [min(wp_tank(:,1)) max(wp_tank(:,1)) ...
                 min(wp_tank(:,2)) max(wp_tank(:,2))] + 60 * [-1 1 -1 1];
% Picked wall/edge points are constrained too, so extend along-wind to them
if ~isempty(picks)
    xy_p = surf_px2world(cal, vertcat(picks.uv));
    xy_p = xy_p(all(isfinite(xy_p), 2) & abs(xy_p(:,1)) < 2000, :);
    if ~isempty(xy_p)
        cal.bbox_mm(1) = min(cal.bbox_mm(1), min(xy_p(:,1)));
        cal.bbox_mm(2) = max(cal.bbox_mm(2), max(xy_p(:,1)));
    end
end

% mm/px at the board centre and at top/bottom of the image (cf. the
% 2.04 / 0.85 mm/px estimated from the wall fits in calib_target_params.m)
cal.mm_per_px_center = mm_per_px(cal, surf_world2px(cal, [opts.along_offset_mm opts.cross_offset_mm]));
cal.mm_per_px_top    = mm_per_px(cal, [cal0.center(1) 1]);
cal.mm_per_px_bottom = mm_per_px(cal, [cal0.center(1) image_size(1)]);

% Distortion must stay monotonic inside the image or the inverse map folds
r_max = max(sqrt(sum(((([1 1; image_size(2) image_size(1)]) - cal.center) / cal.focal_px).^2, 2)));
r = linspace(0, 1.2 * r_max, 200);
slope = ones(size(r));
for i = 1:numel(cal.k), slope = slope + (2*i + 1) * cal.k(i) * r.^(2*i); end
if any(slope <= 0)
    warning(['Fitted radial model is non-monotonic within the image ' ...
        '(k = [%s]); the map is unreliable away from the board.'], sprintf('%.4f ', cal.k));
end

fprintf('\n=== Surface calibration (%s) ===\n', snap_name);
fprintf('Reprojection RMS, homography only : %.3f px\n', rms_h);
fprintf('Reprojection RMS, + radial k1..k%d : %.3f px\n', numel(cal.k), cal.rms_px);
if ~isempty(picks)
    fprintf('  (corners-only fit was %.3f px; lines trade a little corner fit for the periphery)\n', ...
        rms_corners_only);
    for i = 1:numel(picks)
        d = line_dist(cal, picks(i));
        fprintf('  %-22s %2d pts  rms %.2f px  max %.2f px\n', picks(i).name, ...
            numel(d), sqrt(mean(d.^2)), max(d));
    end
end
fprintf('k = [%s], centre = (%.1f, %.1f) px, f = %.0f px\n', ...
    sprintf('%.5f ', cal.k), cal.center, cal.focal_px);
fprintf('Board -> tank axes: [%d %d; %d %d]\n', S_best');
fprintf('mm/px  top %.2f | board centre %.2f | bottom %.2f\n', ...
    cal.mm_per_px_top, cal.mm_per_px_center, cal.mm_per_px_bottom);
fprintf('Trusted footprint (corners + 60 mm): x (along) %.0f..%.0f mm, y (cross) %.0f..%.0f mm\n', ...
    cal.bbox_mm);

%% ---- Check figure ----
if opts.show
    figure('Name', ['Surface calibration ' snap_name], 'Position', [100 100 1300 550]);

    subplot(1,2,1);
    imagesc(frame); colormap(gca, gray); axis image; hold on;
    plot(pts(:,1), pts(:,2), 'g+', 'MarkerSize', 7);
    pr = surf_world2px(cal, wp_tank);
    quiver(pts(:,1), pts(:,2), 20*(pr(:,1)-pts(:,1)), 20*(pr(:,2)-pts(:,2)), 0, 'r');
    overlay_grid(cal, 50, opts.wall_y_mm);
    for i = 1:numel(picks)
        c = line_curve(cal, picks(i));
        plot(c(:,1), c(:,2), 'm-', 'LineWidth', 1);
        plot(picks(i).uv(:,1), picks(i).uv(:,2), 'mo', 'MarkerSize', 5);
    end
    title(sprintf(['Corners (+), residual x20 (red), 50 mm grid, walls (y), ' ...
        'picked lines (m)  rms %.2f px'], cal.rms_px));
    xlabel('u (px)'); ylabel('v (px)');

    subplot(1,2,2);
    [R, xg, yg] = rectify_preview(cal, frame, 1, opts.wall_y_mm + 50);
    imagesc(yg, xg, R); colormap(gca, gray); axis image; hold on;
    plot(wp_tank(:,2), wp_tank(:,1), 'g+', 'MarkerSize', 7);
    xline(opts.wall_y_mm * [-1 1], 'y-', 'LineWidth', 1);
    xlabel('y cross-wind (mm)'); ylabel('x along-wind (mm)');
    title('Rectified onto the surface plane (1 mm grid)');
end

if ~isempty(opts.out_path)
    save(opts.out_path, 'cal');
    fprintf('Saved %s\n', opts.out_path);
end
end


%% =====================================================================
function [pts, bs] = detect_board(frame, expected, flat)
% IR contrast of the anodised/bare-Al board is modest; try a few contrast
% normalisations before giving up.
%
% The cooled detector sees its own reflection (narcissus) as a dark disc +
% ring at the image centre, and vignetting darkens the periphery. A global
% stretch leaves the checker contrast inside the disc too low and the
% detector stops at its edge. Local mean/std normalisation over ~1 square
% removes both, since they vary slowly compared with the 25-45 px squares.
% HighDistortion mode handles the strong barrel bending near the bottom.
if nargin < 3 || isempty(flat), flat = 0; end
f = double(frame) - double(flat);
best_n = 0;
for sigma = [20 30 15 45]
    u8 = local_normalize(f, sigma);
    for high_dist = [false true]
        [pts, bs] = detectCheckerboardPoints(u8, ...
            'PartialDetections', false, 'HighDistortion', high_dist);
        ok = ~isempty(pts) && ~any(isnan(pts(:)));
        if ok && isequal(sort(bs), sort(expected))
            fprintf('Board found (local-norm sigma %d px, HighDistortion %d)\n', ...
                sigma, high_dist);
            return
        end
        if ok && size(pts,1) > best_n
            best_n = size(pts,1); best = {u8, pts, bs, sigma, high_dist};
        end
    end
end
figure('Name', 'Detection failed');
if best_n > 0
    imagesc(best{1}); colormap(gray); axis image; hold on;
    plot(best{2}(:,1), best{2}(:,2), 'r+');
    title(sprintf('Best: %s squares (sigma %d, HighDistortion %d)', ...
        mat2str(best{3}), best{4}, best{5}));
    bs = best{3};
else
    imagesc(local_normalize(f, 20)); colormap(gray); axis image;
end
error('calibrate_surface_single:detection', ...
    ['Expected a %dx%d-square board; best detection %s. See figure. ' ...
     'Try opts.flat (a board-free frame, e.g. mean of the recording) to ' ...
     'subtract the narcissus pattern.'], expected(1), expected(2), mat2str(bs));
end

function u8 = local_normalize(f, sigma)
mu = imgaussfilt(f, sigma);
sd = sqrt(imgaussfilt((f - mu).^2, sigma));
z  = (f - mu) ./ max(sd, eps);
u8 = uint8(255 * (min(max(z, -2.5), 2.5) + 2.5) / 5);
end

function lines = define_lines(wp_tank, sq, opts)
% Lines of known tank position. axis = 1: x = value (spans y);
% axis = 2: y = value (spans x). Pattern edges lie one square beyond the
% outermost inner corners. +y is image right, so the left wall is -y.
xe = [min(wp_tank(:,1)) max(wp_tank(:,1))] + sq * [-1 1];
ye = [min(wp_tank(:,2)) max(wp_tank(:,2))] + sq * [-1 1];
if opts.along_sign > 0, x_top = xe(1); x_bot = xe(2); else, x_top = xe(2); x_bot = xe(1); end
w  = opts.wall_y_mm;
xs = xe + 150 * [-1 1];             % walls extend past the board
ys = ye + 30 * [-1 1];
lines = struct( ...
    'name',  {'left wall', 'right wall', 'pattern top edge', ...
              'pattern bottom edge', 'pattern left edge', 'pattern right edge'}, ...
    'axis',  {2, 2, 1, 1, 2, 2}, ...
    'value', {-w, w, x_top, x_bot, ye(1), ye(2)}, ...
    'span',  {xs, xs, ys, ys, xe, xe}, ...
    'image', {'wall', 'wall', 'snap', 'snap', 'snap', 'snap'}, ...
    'hint',  {'the WATERLINE on the left wall (not its top edge)', ...
              'the WATERLINE on the right wall (not its top edge)', ...
              'the outer edge of the TOP row of squares (incl. its corners)', ...
              'the outer edge of the BOTTOM row of squares (incl. its corners)', ...
              'the outer edge of the LEFT column of squares', ...
              'the outer edge of the RIGHT column of squares'});
end

function picks = get_line_picks(lines, frame, opts, cal)
% Load saved clicks, or collect them with ginput. Only names + pixel
% positions are saved; line geometry is rebuilt from the board each run.
%
% opts.repick: false = reuse saved clicks (click only if none saved),
%              true  = redo all lines,
%              {'pattern bottom edge', 'left wall'} = redo just those;
%              the other lines keep their saved clicks.
picks = [];
saved = struct('name', {}, 'uv', {});
if isfile(opts.picks_path)
    S = load(opts.picks_path, 'saved');
    saved = S.saved;
    fprintf('Loaded line picks from %s\n', opts.picks_path);
end
if iscell(opts.repick) || ischar(opts.repick)
    todo = find(ismember({lines.name}, cellstr(opts.repick)));
    bad = setdiff(cellstr(opts.repick), {lines.name});
    if ~isempty(bad)
        error('Unknown line name(s): %s. Valid: %s', strjoin(bad, ', '), ...
            strjoin({lines.name}, ', '));
    end
elseif opts.repick || isempty(saved)
    todo = 1:numel(lines);
else
    todo = [];
end

if ~isempty(todo)
    % Show the latest full fit (corners + lines) when one is saved, so the
    % dashed line is the estimate you are correcting; else corner-only.
    guess = cal;
    if ~isempty(opts.out_path) && isfile(opts.out_path)
        P = load(opts.out_path, 'cal');
        % Same snapshot -> same corners -> same tank axes, so the saved
        % fit is directly comparable
        if isfield(P.cal, 'snap_path') && strcmp(P.cal.snap_path, opts.snap_path) ...
                && isequal(P.cal.image_size, cal.image_size)
            guess = P.cal;
        end
    end
    fig = figure('Name', 'Pick lines', 'Position', [100 50 1000 820]);
    for n = 1:numel(todo)
        i = todo(n);
        if strcmp(lines(i).image, 'wall'), img = opts.wall_image; else, img = frame; end
        clf(fig);
        imagesc(img); colormap(gray); axis image; hold on;
        s = sort(img(isfinite(img)));
        caxis(double(s(round([0.01 0.99] * numel(s))))');
        c = line_curve(guess, lines(i));
        plot(c(:,1), c(:,2), 'm--', 'LineWidth', 1);
        old = find(strcmp({saved.name}, lines(i).name));
        for j = old
            plot(saved(j).uv(:,1), saved(j).uv(:,2), 'ro', 'MarkerSize', 5);
        end
        auto = [];
        if strcmp(lines(i).image, 'wall') && opts.auto_walls
            auto = auto_wall_points(img, line_curve(cal, lines(i)), opts.wall_band_px);
            plot(auto(:,1), auto(:,2), 'c.', 'MarkerSize', 10);
            hint2 = sprintf(['Cyan = %d auto waterline points. Enter = accept them; ' ...
                'or click your own (replaces)'], size(auto,1));
        elseif ~isempty(old)
            hint2 = 'Red = saved clicks. Enter = KEEP them; or click new ones (replaces)';
        else
            hint2 = 'Enter = done (Enter with no clicks = skip this line)';
        end
        title({sprintf('[%d/%d] Click points along %s', n, numel(todo), lines(i).hint), ...
            ['Magenta dashed = current estimate. Click only where it is off if you like. ' ...
             hint2]}, 'Interpreter', 'none');
        [u, v] = ginput();
        if isempty(u) && ~isempty(auto)
            u = auto(:,1); v = auto(:,2);
        end
        if ~isempty(u)
            saved(old) = [];
            saved(end+1) = struct('name', lines(i).name, 'uv', [u v]); %#ok<AGROW>
        end
    end
    close(fig);
    save(opts.picks_path, 'saved');
    fprintf('Saved line picks to %s\n', opts.picks_path);
end
for j = 1:numel(saved)
    k = find(strcmp({lines.name}, saved(j).name), 1);
    if isempty(k) || size(saved(j).uv, 1) < 2, continue; end
    p = lines(k);
    p.uv = saved(j).uv;
    if isempty(picks), picks = p; else, picks(end+1) = p; end %#ok<AGROW>
end
end

function uv = auto_wall_points(img, c, band)
% Threshold-free waterline finder: along the predicted wall curve c (from
% the corner-only fit), take the strongest cross-wise temperature step
% within +/-band px in each row. No wallTemp/tolerance to tune, unlike the
% commented-out automatic path in distortion_correction.m. Outliers are
% dropped against a quadratic u(v), the same form as the old wall polyfits.
[H, W] = size(img);
[gx, ~] = gradient(imgaussfilt(double(img), 2));
in = c(:,1) >= 1 & c(:,1) <= W & c(:,2) >= 1 & c(:,2) <= H;
c = c(in, :);
uv = zeros(0, 2); strength = zeros(0, 1);
if size(c, 1) < 2, return; end
for v = ceil(min(c(:,2))) : 4 : floor(max(c(:,2)))
    [~, i] = min(abs(c(:,2) - v));
    us = max(1, round(c(i,1) - band)) : min(W, round(c(i,1) + band));
    [m, j] = max(abs(gx(v, us)));
    uv(end+1, :) = [us(j) v];  %#ok<AGROW>
    strength(end+1, 1) = m;    %#ok<AGROW>
end
keep = strength > 0.3 * median(strength);
uv = uv(keep, :);
if size(uv, 1) >= 5
    for pass = 1:2
        p = polyfit(uv(:,2), uv(:,1), 2);
        r = uv(:,1) - polyval(p, uv(:,2));
        uv = uv(abs(r) <= max(3 * 1.4826 * median(abs(r - median(r))), 1.5), :);
    end
end
end

function c = joint_cost(theta, pts, wp, cal, fit_center, picks, w)
cal.k = theta(1:numel(cal.k));
if fit_center, cal.center = theta(end-1:end); end
cal = refit_T(cal, pts, wp);
res = surf_world2px(cal, wp) - pts;
d = line_dist(cal, picks);
c = sqrt(mean([sum(res.^2, 2); (w * d).^2]));
end

function d = line_dist(cal, picks)
% Pixel distance from each picked point to its tank line projected into
% the raw image. Uses only the forward map (always well defined), so it
% stays robust where the inverse iteration would struggle.
d = zeros(0, 1);
for i = 1:numel(picks)
    d = [d; point_to_polyline(picks(i).uv, line_curve(cal, picks(i)))]; %#ok<AGROW>
end
end

function c = line_curve(cal, L)
s = linspace(L.span(1), L.span(2), 400)';
if L.axis == 1
    xy = [L.value * ones(size(s)), s];
else
    xy = [s, L.value * ones(size(s))];
end
c = surf_world2px(cal, xy);
end

function d = point_to_polyline(P, C)
A = C(1:end-1, :);
AB = C(2:end, :) - A;
L2 = max(sum(AB.^2, 2), eps);
d = zeros(size(P, 1), 1);
for i = 1:size(P, 1)
    t = min(max(sum((P(i,:) - A) .* AB, 2) ./ L2, 0), 1);
    d(i) = sqrt(min(sum((A + t .* AB - P(i,:)).^2, 2)));
end
end

function rms = reproj_rms(theta, pts, wp, cal, fit_center)
cal.k = theta(1:numel(cal.k));
if fit_center, cal.center = theta(end-1:end); end
cal = refit_T(cal, pts, wp);
res = surf_world2px(cal, wp) - pts;
rms = sqrt(mean(sum(res.^2, 2)));
end

function cal = refit_T(cal, pts, wp)
% Undistort the detected corners with the current k, then fit H linearly.
cal.T = eye(3);
[~, pu] = surf_px2world(cal, pts);
tf = fitgeotrans(wp, pu, 'projective');
cal.T = tf.T;
end

function s = mm_per_px(cal, uv)
xy = surf_px2world(cal, [uv; uv + [1 0]; uv + [0 1]]);
s = sqrt(abs(det([xy(2,:) - xy(1,:); xy(3,:) - xy(1,:)])));   % sqrt of px area
end

function [uu, vv] = border_pixels(sz)
H = sz(1); W = sz(2);
u = (1:4:W)'; v = (1:4:H)';
uu = [u; W*ones(size(v)); flipud(u); ones(size(v))];
vv = [ones(size(u)); v; H*ones(size(u)); flipud(v)];
end

function overlay_grid(cal, spacing, wall_y)
bb = cal.bbox_mm;
gx = floor(bb(1)/spacing)*spacing : spacing : ceil(bb(2)/spacing)*spacing;
gy = floor(bb(3)/spacing)*spacing : spacing : ceil(bb(4)/spacing)*spacing;
xs = linspace(bb(1), bb(2), 200)'; ys = linspace(bb(3), bb(4), 200)';
for y = gy
    p = surf_world2px(cal, [xs, y*ones(size(xs))]);
    plot(p(:,1), p(:,2), 'c-', 'LineWidth', 0.4);
end
for x = gx
    p = surf_world2px(cal, [x*ones(size(ys)), ys]);
    plot(p(:,1), p(:,2), 'c-', 'LineWidth', 0.4);
end
for y = wall_y * [-1 1]
    p = surf_world2px(cal, [xs, y*ones(size(xs))]);
    plot(p(:,1), p(:,2), 'y-', 'LineWidth', 1.5);
end
end

function [R, xg, yg] = rectify_preview(cal, frame, dx, y_lim)
bb = cal.bbox_mm;
xg = ceil(bb(1)) : dx : floor(bb(2));
yg = max(ceil(bb(3)), -y_lim) : dx : min(floor(bb(4)), y_lim);
[Y, X] = meshgrid(yg, xg);
p = surf_world2px(cal, [X(:) Y(:)]);
R = interp2(double(frame), reshape(p(:,1), size(X)), reshape(p(:,2), size(X)), ...
    'linear', NaN);
end
