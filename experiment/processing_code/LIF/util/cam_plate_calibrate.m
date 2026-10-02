function cal = cam_plate_calibrate(P, opts)
%CAM_PLATE_CALIBRATE Pixel <-> plane map from ONE checkerboard-plate image.
%
%   cal = CAM_PLATE_CALIBRATE(P, opts)
%
% Same model as the IR calibration (IR Camera/calibration/
% calibrate_surface_single.m), so cal works with its surf_world2px /
% surf_px2world:
%
%   p_u = H * [x; y; 1]                                (mm -> undistorted px)
%   p_d = c + (p_u - c) .* (1 + k1 r^2 + k2 r^4 + ...),  r = |p_u - c| / f
%
% fit by minimising corner reprojection error in raw pixels (H solved
% linearly inside, fminsearch over k).  General for the CoreView cameras:
% no wall picks, any 8/16-bit image.
%
% World axes (mm) on the measurement plane, origin at the board centre:
%   x = the board axis that is more horizontal in the image, + to the right
%   y = the other axis, + toward image UP (opts.y_up = true, side views)
%       or image DOWN (opts.y_up = false)
%
% Plate offset: the plate stood opts.offset_mm TOWARD the camera from the
% measurement plane (n_medium: seen through water).  The plane is then
% farther away, so world coordinates are scaled by (D+d)/D about the point
% on the optical axis, D = f_lens_px * mm/px at the board.  Exact for a
% face-on plate, approximate for oblique views.
%
% P    - plate image (2-D; colour images: convert to gray first)
% opts fields (defaults):
%   square_mm   []      measured square size; [] = design of the detected board
%   n_radial    0       radial coefficients. 0 = homography only: the CoreView
%                       lenses are not fisheye (Water_SURF: a uniform scale
%                       alone fits the corners to 0.42 px).  The k1..k2 fit is
%                       still reported as a check that distortion is negligible.
%   lens_f_mm   35      lens focal length (offset scaling + r normalisation)
%   pixel_um    5.5     pixel pitch
%   offset_mm   0       plate toward the camera from the measurement plane
%   n_medium    1.33
%   y_up        true
%   margin_mm   []      trusted footprint = corners +/- margin ([] = 1.5 squares)
%
% cal fields: T, center, focal_px, k, image_size, image_points,
%   world_points, board, square_mm, rms_px, rms_px_homography_only,
%   bbox_mm, mm_per_px_center, offset_scale

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'IR Camera', 'calibration'));   % surf_world2px/px2world, calib_target_params

    def = struct('square_mm', [], 'n_radial', 0, 'lens_f_mm', 35, 'pixel_um', 5.5, ...
                 'offset_mm', 0, 'n_medium', 1.33, 'y_up', true, 'margin_mm', []);
    if nargin < 2, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end

    %% Detect corners
    P   = double(P);
    lim = prctile(P(:), [0.5 99.5]);
    P8  = uint8(255 * min(max((P - lim(1)) / diff(lim), 0), 1));
    [pts, bs] = detectCheckerboardPoints(P8, 'PartialDetections', false);
    if isempty(pts) || any(isnan(pts(:)))
        error('cam_plate_calibrate:detect', 'Checkerboard not detected.');
    end
    board = [];
    for nm = {'fine', 'coarse'}
        b = calib_target_params(nm{1});
        if isequal(sort(bs), sort([b.n_rows b.n_cols])), board = b; end
    end
    if isempty(board)
        error('cam_plate_calibrate:board', 'Detected %dx%d squares: neither fine nor coarse.', bs(1), bs(2));
    end
    if isempty(opts.square_mm), opts.square_mm = board.square_mm; end
    if isempty(opts.margin_mm), opts.margin_mm = 1.5 * opts.square_mm; end
    wp = generateCheckerboardPoints(bs, opts.square_mm);
    wp = wp - mean(wp, 1);

    %% Fit k with H solved linearly inside
    sz = size(P);
    cal.center     = ([sz(2) sz(1)] + 1) / 2;
    cal.focal_px   = opts.lens_f_mm / (opts.pixel_um / 1000);
    cal.k          = zeros(1, opts.n_radial);
    cal.T          = eye(3);
    cal.image_size = sz;

    rms_h = reproj_rms(zeros(1, 0), pts, wp, cal);
    fitk  = @(n) fminsearch(@(th) reproj_rms(th, pts, wp, cal), zeros(1, n), ...
                 optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 5000, 'MaxIter', 5000));
    k2 = fitk(2);                                      % check only
    rms_k2 = reproj_rms(k2, pts, wp, cal);
    if opts.n_radial == 0, k = zeros(1, 0); elseif opts.n_radial == 2, k = k2; else, k = fitk(opts.n_radial); end
    cal.k = k;

    %% Board axes -> world axes from the image orientation
    cal = refit_T(cal, pts, wp);
    p0 = surf_world2px(cal, [0 0]);
    da = surf_world2px(cal, [1 0]) - p0;
    db = surf_world2px(cal, [0 1]) - p0;
    ysgn = -1; if ~opts.y_up, ysgn = 1; end            % image v grows downward
    if abs(da(1)) >= abs(db(1))                        % board a is the more horizontal
        S = [sign(da(1)) 0; 0 ysgn*sign(db(2))];
    else
        S = [0 sign(db(1)); ysgn*sign(da(2)) 0];
    end
    wp = (S * wp.').';

    %% Plate -> measurement-plane offset (scale about the optical axis)
    cal = refit_T(cal, pts, wp);
    xc  = surf_px2world(cal, cal.center);              % plate point on the optical axis
    mmpp_plate = mm_per_px(cal, cal.center);
    D = cal.focal_px * mmpp_plate;
    d = opts.offset_mm / opts.n_medium;
    s = (D + d) / D;
    wp = xc + (wp - xc) * s;
    cal = refit_T(cal, pts, wp);

    res = surf_world2px(cal, wp) - pts;
    cal.rms_px                 = sqrt(mean(sum(res.^2, 2)));
    cal.rms_px_homography_only = rms_h;
    cal.image_points  = pts;
    cal.world_points  = wp;
    cal.board         = board.name;
    cal.board_size    = bs;
    cal.square_mm     = opts.square_mm;
    cal.offset_scale  = s;
    cal.camera_dist_mm = D;
    cal.bbox_mm = [min(wp(:,1)) max(wp(:,1)) min(wp(:,2)) max(wp(:,2))] + opts.margin_mm * [-1 1 -1 1];
    cal.mm_per_px_center = mm_per_px(cal, surf_world2px(cal, [0 0]));
    cal.opts = opts;

    fprintf('Plate: %d corners, ''%s'' board, square %.3f mm\n', size(pts,1), board.name, opts.square_mm);
    fprintf('  reprojection RMS: homography only %.3f px; with radial k1,k2 %.3f px (k = %s) -> using n_radial = %d\n', ...
            rms_h, rms_k2, mat2str(k2, 3), opts.n_radial);
    fprintf('  camera-plate ~%.0f mm; plane offset %.1f mm (apparent %.1f) -> x%.4f\n', ...
            D, opts.offset_mm, d, s);
    fprintf('  %.4f mm/px at the board centre; trusted footprint x %.0f..%.0f, y %.0f..%.0f mm\n', ...
            cal.mm_per_px_center, cal.bbox_mm);
end

function rms = reproj_rms(k, pts, wp, cal)
    cal.k = k;
    cal = refit_T(cal, pts, wp);
    res = surf_world2px(cal, wp) - pts;
    rms = sqrt(mean(sum(res.^2, 2)));
end

function cal = refit_T(cal, pts, wp)
    % Undistort the detected corners with the current k, then fit H linearly
    cal.T = eye(3);
    [~, pu] = surf_px2world(cal, pts);
    tf = fitgeotrans(wp, pu, 'projective');
    cal.T = tf.T;
end

function s = mm_per_px(cal, uv)
    xy = surf_px2world(cal, [uv; uv + [1 0]; uv + [0 1]]);
    s = sqrt(abs(det([xy(2,:) - xy(1,:); xy(3,:) - xy(1,:)])));   % sqrt of px area
end
