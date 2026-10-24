function cal = cisg_plate_calibrate(cfg)
%CISG_PLATE_CALIBRATE Pixel <-> tank-plane map for the nadir CISG camera from
% ONE view of the big checkerboard plate, in the SAME tank frame as the IR cal.
%
%   cal = CISG_PLATE_CALIBRATE            % defaults (CoreView 97 plate frame)
%   cal = CISG_PLATE_CALIBRATE(cfg)
%
% The model is the IR one (IR Camera/calibration/calibrate_surface_single.m),
% so cal works with surf_world2px / surf_px2world:
%
%   p_u = H * [x; y; 1]                                (mm -> undistorted px)
%   p_d = c + (p_u - c) .* (1 + k1 r^2 + k2 r^4 + ...),  r = |p_u - c| / f
%
% fit by minimising corner reprojection error in raw pixels (H solved
% linearly inside, fminsearch over k).
%
% Frame.  The plate lay AT the water surface for this frame and for IR
% Snap-000096 (taken at the same moment), so both cals describe the surface
% plane.  The lattice is identified in the IR cal's PLATE frame (origin at the
% plate centre; IR cal: tank = board_to_tank * board), then rotated by
% theta_deg about the plate centre into the CAMERA-aligned tank frame stored in
% cal (mm): x = image right = downwind, y = image up = cross-wind.  The camera
% is aligned with the wind, the plate may be slightly rotated.
% cal.plate_to_tank (2x2, columns) maps plate-frame mm to tank mm.
%
% COORDINATE CONVENTION (tank frame, cm in the slope outputs).
%   x: downwind = +image right (camera aligned with the wind); x = 0 on the laser-sheet line, which is the plate
%      centre line along x (cfg.laser_x_px picks it; the IR cal's board centre is the true plate centre, so the
%      lattice is placed with that corner line at x = 0).  The IR along-wind axis (image top -> bottom) must map
%      to +x: check with cisg_check_colocation (IR +x direction vs CISG +x) after rebuilding with use_ir.
%   y: +image up (cross-wind, handedness unchanged); the internal frame has y = 0 at a lattice line guessed
%      from the image centre.  The channel centre y0 (wall edges of the run's reference, cisg_channel_center.m)
%      is applied only at OUTPUT (X_cm/Y_cm, plots: y_out = y - y0), so IR <-> CISG mapping stays in the plate frame.
%
% The CISG sees only ~7 x 6 of the 10 x 9 inner corners (0.07 mm/px, 292 x
% 219 mm), so detectCheckerboardPoints returns a SUB-BLOCK and its corner
% ordering is that of the block, not of the board.  The block is therefore
% placed on the board (rather than assuming the IR corner ordering) by:
%   1. axes:    tank x is the image-horizontal axis (thick-border axis, see
%               make_calibration_target_svg.m) and the view is not mirrored
%               (same handedness as IR) -> two options, 180 deg apart;
%   2. offset:  the plate centre is taken to be the image centre (nearest whole
%               square; the lattice position of the image centre is checked to
%               be a corner line along x and a square centre along y);
%   3. sign:    the 180-deg pair is split by square parity.  Ablated (bare Al)
%               squares are bright in IR; in Snap-000096 they are the squares
%               with odd (m + n), m = round(x/40 - 1/2), n = round(y/40).  They
%               are assumed bright in the CISG image too (cfg.ablated_bright).
%
% cfg fields (defaults):
%   frame_path    CoreView_97 ..._001.raw   plate frame (colour, gray = rgb2gray)
%   ir_cal_path   IR Camera/calibration/surface_cal_Snap-000096.mat
%   out_path      calibration/cisg_cal_CoreView_97.mat ('' = no save)
%   board         'coarse'  calib_target_params name
%   square_mm     []        measured corner pitch; [] = design
%   scales        [0.5 0.25 0.125]  detection image scales tried
%   sigma_px      [40 20 80]        local-normalisation sigma at full res (px)
%   n_radial      1         radial coefficients (k1..kn); the fit with k1,k2 is
%                           still reported
%   lens_f_mm     35        lens focal length (only normalises r)
%   pixel_um      5.5       pixel pitch
%   ablated_bright true     SVG-black (engraved, anodising removed) squares are bright in the CISG image
%   laser_x_px    NaN       x = 0 on the laser-sheet line: NaN = not set (x origin = lattice node nearest the image centre,
%                           a guess that can be off by whole squares), [] = pick it interactively on the plate frame,
%                           or an image column (px); the lattice is re-placed so the nearest corner line is x = 0
%   use_ir        true      false = camera-only frame (new plate frames without an IR snap): tank x = image right,
%                           y = image up, origin = plate lattice node nearest the image centre; no IR cal, no parity
%                           or thick-axis checks (those need the IR snapshot)
%   ablated_parity 1        parity of (m + n) of the ablated squares in the IR tank frame
%   show          false
%
% cal fields: T, center, focal_px, k, image_size, image_points, world_points,
%   board_to_tank, rms_px, mm_per_px_center, plate_centre_px,
%   plate_centre_offset_px, bbox_mm, checks, frame_path, ir_cal_path, cfg, ...

    here = fileparts(mfilename('fullpath'));
    sg   = fileparts(here);
    addpath(sg, fullfile(fileparts(sg), 'IR Camera', 'calibration'));   % cisg_load_coreview, surf_*, calib_target_params

    def = struct('frame_path', fullfile('D:\HLAB_2026\Core', 'CoreView_97', 'Flare 12M125 CCL C1112A00', ...
                                        'CoreView_97_Flare 12M125 CCL C1112A00_001.raw'), ...
                 'ir_cal_path', fullfile(fileparts(sg), 'IR Camera', 'calibration', 'surface_cal_Snap-000096.mat'), ...
                 'out_path', fullfile(here, 'cisg_cal_CoreView_97.mat'), ...
                 'board', 'coarse', 'square_mm', [], 'scales', [0.5 0.25 0.125], 'sigma_px', [40 20 80], ...
                 'n_radial', 1, 'lens_f_mm', 35, 'pixel_um', 5.5, ...
                 'ablated_bright', true, 'ablated_parity', 1, 'use_ir', true, 'frame_size', [3072 4100], 'laser_x_px', NaN, 'orientation_checked', false, 'orientation_rule', '', 'show', false);
    if nargin < 1, cfg = struct(); end
    for f = fieldnames(def).'
        if ~isfield(cfg, f{1}), cfg.(f{1}) = def.(f{1}); end
    end
    board = calib_target_params(cfg.board);
    if isempty(cfg.square_mm), cfg.square_mm = board.square_mm; end
    sq = cfg.square_mm;

    %% IR cal: board -> tank axes
    if cfg.use_ir && ~isfile(cfg.ir_cal_path)
        error('cisg_plate_calibrate:ir_cal', 'IR calibration not found: %s', cfg.ir_cal_path);
    end
    cal_ir = [];  S_ir = [];  off_ir = [0 0];
    if cfg.use_ir
    S = load(cfg.ir_cal_path, 'cal');
    cal_ir = S.cal;
    S_ir = cal_ir.board_to_tank;                       % tank = S_ir * board (+ offsets)
    off_ir = [cal_ir.opts.along_offset_mm cal_ir.opts.cross_offset_mm];
    if any(off_ir ~= 0) || cal_ir.opts.along_sign ~= 1
        warning('cisg_plate_calibrate:ir_offsets', ...
            'IR cal uses offsets / along_sign (%s, %d); the parity rule below was measured for 0 / +1.', ...
            mat2str(off_ir), cal_ir.opts.along_sign);
    end
    end

    %% Detect the corner block in the CISG frame
    I = double(rgb2gray(cisg_load_coreview(cfg.frame_path, cfg.frame_size)));   % frame_size: headerless files
    sz = size(I);
    [pts, bs] = detect_block(I, cfg);
    nr = bs(1) - 1;  nc = bs(2) - 1;                   % inner corners, y fastest (detector order)
    fprintf('CISG plate frame: %d x %d px, block %d x %d corners (%d)\n', sz(2), sz(1), nc, nr, size(pts,1));

    %% Place the block on the board (axes, offset, sign)
    [ix, iy] = meshgrid(0:nc-1, 0:nr-1);
    L  = [ix(:) iy(:)];                                % lattice units, integers at corners
    c0 = ([sz(2) sz(1)] + 1) / 2;
    tf = fitgeotrans(L, pts, 'projective');
    cl = transformPointsInverse(tf, c0);               % image centre in lattice units
    Jp = transformPointsForward(tf, cl + [0.5 0]) - transformPointsForward(tf, cl - [0.5 0]);
    Jq = transformPointsForward(tf, cl + [0 0.5]) - transformPointsForward(tf, cl - [0 0.5]);
    J  = [Jp' Jq'];                                    % d(u,v)/d(lattice p,q)

    % x horizontal, not mirrored (IR: x ~ v, y ~ u has det < 0) -> two options
    opt = struct('R', {}, 'sh', {}, 'wp', {}, 'centre_mm', {}, 'contrast', {});
    for P = {[1 0; 0 1], [0 1; 1 0]}
        for sgn = [1 1; 1 -1; -1 1; -1 -1]'
            R  = diag(sgn) * P{1};                     % tank(squares) = R * lattice
            Jt = J * R';                               % columns d(u,v)/dx, d(u,v)/dy
            if abs(Jt(1,1)) > abs(Jt(2,1)) && det(Jt) < 0
                q  = (R * (L - mean(L, 1))')';
                cq = (R * (cl - mean(L, 1))')';
                sx = round(q(1,1)) - q(1,1);           % x on whole squares (corner lines)
                sy = 0.5 - q(1,2) - floor(0.5 - q(1,2));   % y on half squares (square centres)
                cp = cq + [sx sy];
                sh = [sx sy] - round(cp);              % plate centre -> nearest image centre
                opt(end+1).R = R;                      %#ok<AGROW>
                opt(end).sh = sh;
                opt(end).wp = sq * (q + sh);
                opt(end).centre_mm = sq * (cq + sh);   % image centre in tank mm
                opt(end).jxx = Jt(1,1);                % > 0: tank x points image right
            end
        end
    end
    if numel(opt) ~= 2
        error('cisg_plate_calibrate:axes', 'Expected 2 axis options (x horizontal, unmirrored), got %d.', numel(opt));
    end

    % Image-centre lattice position vs the plate-centre expectation:
    % along x a corner line (A), along y a square centre; B is the transpose
    ah = abs(Jp(1)) > abs(Jp(2));                      % lattice p is the horizontal axis
    dev = abs(cl - round(cl));                         % 0 = on a corner line, 0.5 = square centre
    if ~ah, dev = dev([2 1]); end                      % dev(1): horizontal, dev(2): vertical
    misfit_A = sq * (dev(1) + (0.5 - dev(2)));
    misfit_B = sq * ((0.5 - dev(1)) + dev(2));

    % Parity: ablated squares should be brighter (cfg.ablated_bright)
    Is = imgaussfilt(I, 40);                           % ~ 0.07 square
    [cx, cy] = meshgrid(0.5:nc-1.5, 0.5:nr-1.5);       % cell centres, lattice units
    cpx = transformPointsForward(tf, [cx(:) cy(:)]);
    bright = interp2(Is, cpx(:,1), cpx(:,2), 'linear', NaN);
    for i = 1:numel(opt)
        q  = (opt(i).R * ([cx(:) cy(:)] - mean(L, 1))')';
        xyc = sq * (q + opt(i).sh);
        m = round(xyc(:,1) / sq - 0.5);  n = round(xyc(:,2) / sq);
        ab = mod(m + n, 2) == cfg.ablated_parity;
        opt(i).contrast = mean(bright(ab), 'omitnan') - mean(bright(~ab), 'omitnan');
    end
    [~, ibest] = max(([opt.contrast]) * (2 * cfg.ablated_bright - 1));
    if ~cfg.use_ir, ibest = find([opt.jxx] > 0, 1); end   % camera frame: x = image right, no parity needed
    wp = opt(ibest).wp;

    % Laser-sheet line = plate centre line = x = 0: the corner line under the laser is re-placed at x = 0
    % (the guessed placement above can be off by whole 40 mm squares).  cfg.laser_x_px: NaN = not set (keep the
    % guess), [] = pick interactively, number = image column; the pick snaps to the nearest corner line.
    laser = struct('px', NaN, 'snapped_px', NaN, 'x_shift_mm', 0);
    if isempty(cfg.laser_x_px) || ~isnan(cfg.laser_x_px)
        xl = unique(round(wp(:,1) * 1e3) / 1e3);                       % x of the detected vertical corner lines (mm)
        ul = arrayfun(@(x) median(pts(abs(wp(:,1) - x) < 1e-3, 1)), xl);   % their image columns
        u = cfg.laser_x_px;
        if isempty(u), u = pick_laser_column(I, ul); end
        [~, il] = min(abs(ul - u));
        laser = struct('px', u, 'snapped_px', ul(il), 'x_shift_mm', xl(il));
        wp(:,1) = wp(:,1) - xl(il);                                    % x = 0 on the laser line
        fprintf('  laser line at column %.0f px -> corner line at %.0f px, shifts the x origin by %+.0f mm (%+.0f squares)\n', ...
                u, ul(il), xl(il), xl(il) / sq);
    end
    fprintf('  axis options (contrast ablated - other, counts): %+.1f | %+.1f -> option %d\n', ...
            opt(1).contrast, opt(2).contrast, ibest);

    %% Fit k with H solved linearly inside
    cal = struct('center', c0, 'focal_px', cfg.lens_f_mm / (cfg.pixel_um / 1000), ...
                 'k', zeros(1, 0), 'T', eye(3), 'image_size', sz);
    rms_h = reproj_rms(zeros(1, 0), pts, wp, cal);
    fitk  = @(n) fminsearch(@(th) reproj_rms(th, pts, wp, cal), zeros(1, n), ...
                 optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 5000, 'MaxIter', 5000));
    k2 = fitk(2);
    rms_k2 = reproj_rms(k2, pts, wp, cal);
    if cfg.n_radial == 0, k = zeros(1, 0); elseif cfg.n_radial == 2, k = k2; else, k = fitk(cfg.n_radial); end
    cal.k = k;
    cal = refit_T(cal, pts, wp);

    res = surf_world2px(cal, wp) - pts;
    cal.rms_px                 = sqrt(mean(sum(res.^2, 2)));
    cal.rms_px_homography_only = rms_h;
    cal.rms_px_k1k2            = rms_k2;
    cal.image_points  = pts;
    cal.world_points  = wp;
    cal.board         = board.name;
    cal.board_size    = bs;
    cal.board_to_tank = S_ir;
    cal.square_mm     = sq;
    cal.bbox_mm = [min(wp(:,1)) max(wp(:,1)) min(wp(:,2)) max(wp(:,2))] + 1.5 * sq * [-1 1 -1 1];
    cal.plate_centre_px        = surf_world2px(cal, [0 0]);
    cal.image_centre_px        = c0;
    cal.plate_centre_offset_px = cal.plate_centre_px - c0;
    cal.mm_per_px_center       = mm_per_px(cal, cal.plate_centre_px);
    cal.frame_path  = cfg.frame_path;
    cal.ir_cal_path = iff_(cfg.use_ir, cfg.ir_cal_path, '');
    cal.orientation_rule = cfg.orientation_rule;  if isempty(cal.orientation_rule), cal.orientation_rule = 'parity (consistent with wind direction)'; end
    cal.laser_x_px = laser.px;  cal.laser_snapped_px = laser.snapped_px;  cal.laser_x_shift_mm = laser.x_shift_mm;
    cal.cfg = cfg;

    %% Plate frame -> camera-aligned tank frame
    % The camera is aligned with the wind, the plate may be rotated.  Tank x =
    % image +u (right, downwind), tank y = image UP (-v): right-handed with y
    % up, the unmirrored (IR-like) handedness.  theta = angle of plate +x from
    % the image rows, counter-clockwise on the display (toward image up).
    % Rotate about the plate centre: tank = plate_to_tank * plate (columns).
    p0  = cal.plate_centre_px;
    Jp0 = [surf_world2px(cal, [1 0]) - p0; surf_world2px(cal, [0 1]) - p0]';  % d(u,v)/d(plate x,y)
    theta = atan2d(-Jp0(2,1), Jp0(1,1));
    if cfg.use_ir && ~cfg.orientation_checked
        % physical constraint: the camera is aligned with the wind and IR +x (along-wind) must map to CISG image right
        if abs(theta) > 135       % parity chose the wrong 180 deg option: take the other one
            fprintf('  parity overridden by wind-direction constraint (theta %+.1f deg, parity contrast %+.1f counts)\n', theta, opt(ibest).contrast);
            cfg2 = cfg;  cfg2.ablated_parity = 1 - cfg.ablated_parity;  cfg2.orientation_checked = true;
            cfg2.orientation_rule = sprintf('parity overridden by wind-direction constraint (parity gave theta %+.1f deg, contrast %+.1f)', theta, opt(ibest).contrast);
            cal = cisg_plate_calibrate(cfg2);  return;
        elseif abs(theta) > 45
            error('cisg_plate_calibrate:axes', 'theta = %+.1f deg: an axis swap, not a parity issue (IR +x does not map to the CISG image right).', theta);
        end
    end
    Rot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    cal.theta_deg       = theta;
    cal.plate_to_tank   = Rot;
    cal.T               = blkdiag(Rot, 1) * cal.T;     % plate_row = tank_row * Rot
    cal.world_points_plate = wp;
    cal.world_points    = wp * Rot';
    cal.bbox_mm = [min(cal.world_points(:,1)) max(cal.world_points(:,1)) ...
                   min(cal.world_points(:,2)) max(cal.world_points(:,2))] + 1.5 * sq * [-1 1 -1 1];
    cal.plate_centre_px = surf_world2px(cal, [0 0]);   % unchanged (origin fixed)

    %% Axis assertions
    Jc = [surf_world2px(cal, [1 0]) - p0; surf_world2px(cal, [0 1]) - p0]';   % d(u,v)/d(tank x,y)
    chk.x_horizontal = abs(Jc(1,1)) > abs(Jc(2,1));            % (a) tank x = image horizontal
    chk.x_image_right = Jc(1,1) > 0;                           % tank +x = image right
    chk.theta_small  = abs(theta) < 45;                        % plate roughly aligned with the camera
    chk.unmirrored   = det(Jc) < 0;                            % same handedness as IR
    % (b) tank x is the thick-border board axis (IR cal board axes)
    thick_n = board.n_rows; if board.board_x0_mm > board.board_y0_mm, thick_n = board.n_cols; end
    nsq_x = thick_n;  chk.x_is_thick_axis = true;              % (no IR cal: not checked)
    if cfg.use_ir
        jx = find(S_ir(1,:) ~= 0, 1);                          % board column that is tank x
        nsq_x = cal_ir.board_size(3 - jx);                     % squares spanned by it (col1 <-> bs(2))
        chk.x_is_thick_axis = (nsq_x == thick_n);
    end
    chk.centre_misfit_mm = [misfit_A misfit_B];                % x-corner-line hypothesis vs transposed
    chk.centre_hypothesis_ok = misfit_A < misfit_B;
    cal.checks = chk;
    if ~chk.x_horizontal, warning('cisg_plate_calibrate:axis', 'Tank x is not the CISG horizontal axis.'); end
    if ~chk.theta_small
        warning('cisg_plate_calibrate:theta', ['Plate +x is rotated %.1f deg from the camera x axis: the IR cal ' ...
            'along/cross convention disagrees with the camera (not patched).'], theta);
    end
    if ~chk.unmirrored,   warning('cisg_plate_calibrate:axis', 'CISG map is mirrored relative to IR.'); end
    if ~chk.x_is_thick_axis
        warning('cisg_plate_calibrate:axis', 'IR tank x spans %d squares; thick-border axis spans %d.', nsq_x, thick_n);
    end
    if ~chk.centre_hypothesis_ok
        warning('cisg_plate_calibrate:axis', ['Image centre sits nearer a square centre than a corner line along ' ...
            'x (misfit %.1f mm vs %.1f mm transposed); check the block placement.'], misfit_A, misfit_B);
    end

    fprintf('  reprojection RMS: homography only %.3f px; k1,k2 %.3f px (k = %s); using n_radial = %d -> %.3f px\n', ...
            rms_h, rms_k2, mat2str(k2, 3), cfg.n_radial, cal.rms_px);
    fprintf('  %.5f mm/px at the plate centre; plate centre at (%.1f, %.1f) px, offset from image centre (%+.1f, %+.1f) px\n', ...
            cal.mm_per_px_center, cal.plate_centre_px, cal.plate_centre_offset_px);
    fprintf('  plate rotation theta = %+.2f deg (plate +x from image rows, CCW on display); tank +x = image right, +y = image up\n', theta);
    fprintf('  checks: x horizontal %d | unmirrored %d | IR x = thick axis (%d sq) %d | centre misfit %.1f mm (transposed %.1f)\n', ...
            chk.x_horizontal, chk.unmirrored, thick_n, chk.x_is_thick_axis, misfit_A, misfit_B);

    if ~isempty(cfg.out_path)
        save(cfg.out_path, 'cal');
        fprintf('Saved %s\n', cfg.out_path);
    end
end

function [pts, bs] = detect_block(I, cfg)
    % Largest corner block over a few image scales and local normalisations
    % (the IR code does the same: the plate contrast varies across the frame).
    % Squares are ~560 px here, so detection runs on reduced images.
    best_n = 0;  pts = [];  bs = [];
    for s = cfg.scales
        Ir = imresize(I, s, 'box');
        for sg = cfg.sigma_px
            u8 = local_normalize(Ir, sg * s);
            ws = warning('off', 'all');                % even x even blocks warn about orientation
            [p, b] = detectCheckerboardPoints(u8, 'PartialDetections', false);
            warning(ws);
            if ~isempty(p) && ~any(isnan(p(:))) && size(p,1) > best_n
                best_n = size(p,1);
                pts = (p - 0.5) / s + 0.5;             % back to full-res pixel coordinates
                bs = b;
            end
        end
    end
    if best_n < 20
        error('cisg_plate_calibrate:detect', 'Checkerboard block not detected (best %d corners).', best_n);
    end
end

function u8 = local_normalize(f, sigma)
    mu = imgaussfilt(f, sigma);
    sd = sqrt(imgaussfilt((f - mu).^2, sigma));
    z  = (f - mu) ./ max(sd, eps);
    u8 = uint8(255 * (min(max(z, -2.5), 2.5) + 2.5) / 5);
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

function v = iff_(c, a, b)
    if c, v = a; else, v = b; end
end

function u = pick_laser_column(I, ul)
    % show the plate frame with the detected vertical corner lines; click the laser-sheet line
    fig = figure('Name', 'Pick the laser line', 'Position', [100 100 1300 950]);
    imagesc(I(1:4:end, 1:4:end)); colormap(gray); axis image; hold on;
    for k = 1:numel(ul), xline(ul(k) / 4, 'c-'); end
    title({'Click the LASER-SHEET line (the plate centre line, vertical in the image); it snaps to the nearest corner line (cyan)'});
    [u, ~] = ginput(1);  u = u * 4;
    close(fig);
end
