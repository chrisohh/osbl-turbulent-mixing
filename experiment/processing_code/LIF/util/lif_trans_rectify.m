function T = lif_trans_rectify(fplate, opts)
%LIF_TRANS_RECTIFY Flatten the oblique transverse view onto the laser-sheet plane.
%
%   T = LIF_TRANS_RECTIFY(fplate, opts)      build the map from a plate image
%   J = lif_trans_apply(T, I)                resample a frame onto the plane
%
% The bottom (Water_PIV) camera sees the laser sheet ~58 deg off its axis,
% so a uniform scale only holds vertically and the cross-tank axis is
% squeezed to ~half.  The plate imaged IN the sheet (run 103 for 104/105,
% run 121's laser-sheet image for 119/121) gives the full map, with the
% same model as the IR fisheye correction (IR Camera/calibration,
% cam_distortion_correction.m):
%
%   p_u = H * [y; z; 1]                                (plate mm -> undistorted px)
%   p_d = c + (p_u - c) .* (1 + k1 r^2 + k2 r^4),      r = |p_u - c| / f
%
% H solved linearly inside, fminsearch over k, minimising the corner
% reprojection error in raw pixels.  Every frame is then resampled
% (surf_world2px + interp2) on a regular (y, z) grid in the plate plane --
% as if the camera looked straight along the plate normal -- with true cm
% on BOTH axes.  T is a cal struct, so surf_world2px / surf_px2world work
% on it directly.
%
% Unlike cam_plate_calibrate, partial boards are fine (the square size is
% given, not looked up from the fine/coarse design).
%
% fplate - plate image file, or a cell list of plate frames (first usable)
% opts (defaults):
%   square_mm 11, n_radial 0 (homography only: the board covers a quarter
%      of the frame, so radial terms are unconstrained outside it and only
%      gain ~0.01 px on it -- the 60 mm lens is not a fisheye; the k1,k2
%      fit is still printed as a check), lens_f_mm 60,
%   pixel_um 5.5 (only normalise r), dx_mm 0.12 (output grid spacing),
%   cols [] (image columns to use, e.g. [1 3900] to drop a bright edge),
%   outlier_px 3 (corners off the homography by more than max(this,
%      3 x median) are dropped as mis-detections),
%   cache_file '' -- the model (not the large U/V maps) is saved here when
%      the plate is found, and loaded from here when it is not (plate
%      images on another computer, as lif_scale_cached),
%   verbose true
%
% T fields: T, center, focal_px, k, image_size (the cal model);
%   y_mm, z_mm (output grid vectors; z decreasing down the rows, origin at
%   the board centre, z up, y toward image right), U, V (raw-image position
%   of every output pixel, NaN outside the footprint), footprint, dx_mm,
%   rms_px, rms_px_homography_only, board_size, file

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'IR Camera', 'calibration'));   % surf_world2px / surf_px2world

    def = struct('square_mm', 11, 'n_radial', 0, 'lens_f_mm', 60, 'pixel_um', 5.5, ...
                 'dx_mm', 0.12, 'cols', [], 'outlier_px', 3, 'cache_file', '', 'verbose', true);
    if nargin < 2, opts = struct(); end
    for f = fieldnames(def).'
        if ~isfield(opts, f{1}), opts.(f{1}) = def.(f{1}); end
    end
    if ischar(fplate) || isstring(fplate), fplate = cellstr(fplate); end

    %% Detect corners (first usable plate frame)
    pts = [];
    for i = 1:numel(fplate)
        if ~isfile(fplate{i}), continue; end
        P   = lif_load_raw(fplate{i});
        lim = prctile(P(:), [0.5 99.5]);
        P8  = uint8(255 * min(max((P - lim(1)) / diff(lim), 0), 1));
        [p, bs] = detectCheckerboardPoints(P8, 'PartialDetections', false);
        if ~isempty(p) && ~any(isnan(p(:))) && all(bs >= 3)
            pts = p;  board_size = bs;  file = fplate{i};  break
        end
        if opts.verbose, fprintf('  %s: no board -- trying the next plate frame\n', fplate{i}); end
    end
    if isempty(pts)
        if ~isempty(opts.cache_file) && isfile(opts.cache_file)
            C = load(opts.cache_file);
            T = build_grid(C.T, opts);
            if opts.verbose
                fprintf('Plate image not here; flattening map from %s (RMS %.2f px, k = %s)\n', ...
                        opts.cache_file, T.rms_px, mat2str(T.k, 3));
            end
            return
        end
        error('lif_trans_rectify:detect', 'No checkerboard found in the plate frames (and no cache).');
    end
    sz = size(P);
    wp = generateCheckerboardPoints(board_size, opts.square_mm);
    wp = wp - mean(wp, 1);                              % origin at the board centre

    %% Drop mis-detected corners (plate 103: a few ~10 px off, against the
    %  glare at the board edge) -- homography residual > max(3 px, 3 x median)
    n0 = size(pts, 1);
    for it = 1:3
        tf = fitgeotrans(wp, pts, 'projective');
        r  = sqrt(sum((transformPointsForward(tf, wp) - pts).^2, 2));
        bad = r > max(opts.outlier_px, 3 * median(r));
        if ~any(bad), break; end
        pts(bad, :) = [];  wp(bad, :) = [];
    end
    if size(pts, 1) < 0.7 * n0
        warning('lif_trans_rectify:outliers', 'Only %d of %d corners kept -- check the plate detection.', size(pts,1), n0);
    end

    %% Fit k with H solved linearly inside (as cam_plate_calibrate)
    T.center     = ([sz(2) sz(1)] + 1) / 2;
    T.focal_px   = opts.lens_f_mm / (opts.pixel_um / 1000);
    T.k          = zeros(1, 0);
    T.T          = eye(3);
    T.image_size = sz;
    rms_h = reproj_rms(zeros(1, 0), pts, wp, T);
    fitk = @(n) fminsearch(@(th) reproj_rms(th, pts, wp, T), zeros(1, n), ...
                 optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 5000, 'MaxIter', 5000));
    k2 = fitk(2);                                       % check only
    rms_k2 = reproj_rms(k2, pts, wp, T);
    if opts.n_radial == 0, k = zeros(1, 0); elseif opts.n_radial == 2, k = k2; else, k = fitk(opts.n_radial); end
    T.k = k;

    %% Board axes -> (y, z): z = the axis more vertical in the image, + up
    T  = refit_T(T, pts, wp);
    p0 = surf_world2px(T, [0 0]);
    da = surf_world2px(T, [1 0]) - p0;
    db = surf_world2px(T, [0 1]) - p0;
    if abs(da(1)) >= abs(db(1))                         % board a is the more horizontal
        S = [sign(da(1)) 0; 0 -sign(db(2))];
    else
        S = [0 sign(db(1)); -sign(da(2)) 0];
    end
    wp = (S * wp.').';
    T  = refit_T(T, pts, wp);
    res = surf_world2px(T, wp) - pts;

    T.rms_px = sqrt(mean(sum(res.^2, 2)));
    T.rms_px_homography_only = rms_h;
    T.board_size = board_size;  T.square_mm = opts.square_mm;
    T.image_points = pts;  T.world_points = wp;  T.file = file;
    if ~isempty(opts.cache_file)
        if ~isfolder(fileparts(opts.cache_file)), mkdir(fileparts(opts.cache_file)); end
        save(opts.cache_file, 'T');
    end
    T = build_grid(T, opts);

    if opts.verbose
        q = surf_px2world(T, [p0; p0 + [1 0]; p0 + [0 1]]);
        fprintf('Rectify %s\n  %d of %d corners (%dx%d squares of %.1f mm); reprojection RMS: homography only %.2f px; with radial k1,k2 %.2f px (k = %s) -> n_radial = %d\n', ...
                file, size(pts,1), n0, board_size, opts.square_mm, rms_h, rms_k2, mat2str(k2, 3), opts.n_radial);
        fprintf('  raw image at the board centre: %.4f mm/px horizontal, %.4f mm/px vertical (ratio %.2f)\n', ...
                norm(q(2,:) - q(1,:)), norm(q(3,:) - q(1,:)), norm(q(2,:) - q(1,:)) / norm(q(3,:) - q(1,:)));
        fprintf('  flattened grid %d x %d at %.3f mm: y %.0f..%.0f mm, z %.0f..%.0f mm (board centre = 0)\n', ...
                numel(T.z_mm), numel(T.y_mm), opts.dx_mm, T.y_mm([1 end]), T.z_mm([end 1]));
    end
end

function T = build_grid(T, opts)
    % Regular (y, z) grid over the footprint of the used image region on
    % the plane, and the raw-image position of each of its pixels
    sz = T.image_size;
    cols = opts.cols;  if isempty(cols), cols = [1 sz(2)]; end
    n = 60;
    edge = [linspace(cols(1), cols(2), n)' ones(n,1); ...
            cols(2)*ones(n,1) linspace(1, sz(1), n)'; ...
            linspace(cols(2), cols(1), n)' sz(1)*ones(n,1); ...
            cols(1)*ones(n,1) linspace(sz(1), 1, n)'];
    yz = surf_px2world(T, edge);
    y_mm = ceil(min(yz(:,1))) : opts.dx_mm : floor(max(yz(:,1)));
    z_mm = floor(max(yz(:,2))) : -opts.dx_mm : ceil(min(yz(:,2)));   % row 1 = top
    [Y, Z] = meshgrid(y_mm, z_mm);
    uv = surf_world2px(T, [Y(:) Z(:)]);
    U = reshape(uv(:,1), size(Y));  V = reshape(uv(:,2), size(Y));
    out = ~(U >= cols(1) & U <= cols(2) & V >= 1 & V <= sz(1));
    U(out) = NaN;  V(out) = NaN;
    T.y_mm = y_mm;  T.z_mm = z_mm;  T.U = U;  T.V = V;  T.footprint = ~out;
    T.dx_mm = opts.dx_mm;  T.cols = opts.cols;
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
