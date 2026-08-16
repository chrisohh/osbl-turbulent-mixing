function p = calib_target_params(name)
%CALIB_TARGET_PARAMS  Shared geometry for the laser-cut calibration targets.
%
%   p = calib_target_params('coarse')   % IR board
%   p = calib_target_params('fine')     % color / mono board
%
% Single source of truth for square size and corner counts. Called by
% make_calibration_target_svg.m, make_synthetic_calib_data.m, and both
% calibration drivers, so the laser-cut plate and corners.csv cannot drift
% apart.
%
% Square counts are deliberately odd x even. detectCheckerboardPoints
% requires one even and one odd edge to resolve board orientation
% unambiguously; an even x even board is 180-deg ambiguous.
%
% Coordinate convention for the emitted corner list (matches
% generateCheckerboardPoints and detectCheckerboardPoints output ordering):
% row-major with Y varying fastest, X slowest, origin at the board centre.

switch lower(name)
    case 'coarse'
        % IR board. 40 mm squares - sized from REAL wall-detection data
        % (Rec-000014, p_left/p_right polyfits), not a theoretical FOV
        % calc. The camera views the water surface at a steep oblique
        % angle, so resolution varies sharply across the frame:
        %   row   1 (far/top)   : 2.04 mm/px
        %   row 512 (near/bottom): 0.85 mm/px   (2.4x range)
        % A uniform square size has to survive the WORST case (far/top),
        % where the ~20 px/square detectability floor requires ~41 mm
        % squares. 20 mm (the earlier test-coupon-validated size) only
        % works in the near part of the frame and would fail detection
        % near the top.
        %
        % The same wall-curvature data also confirms strong lens
        % distortion independent of the perspective/trapezoid effect: a
        % physically straight wall maps to a straight line under an ideal
        % pinhole projection, so the ~0.0007/-0.0005 quadratic terms in
        % the polyfits are pure lens bending (~130-184 px of curvature
        % across the frame height, ~20-29% of image width) - see
        % calibrate_intrinsic.m's NumRadialDistortionCoefficients=3 and
        % its radial-residual diagnostic.
        p.name        = 'coarse';
        p.square_mm   = 40;
        p.n_cols      = 11;      % squares across  (odd)  -> 440 mm
        p.n_rows      = 10;      % squares down    (even) -> 400 mm
        p.quiet_mm    = 10;      % ablated border ring width - kept narrow
                                  % so width-axis margin to the plate edge
                                  % stays a comfortable ~11 mm (was 1.3 mm
                                  % at quiet_mm=20, uncomfortably close to
                                  % laser cut-tolerance/stock-size error)
        p.plate_w_mm  = 482.6;   % 19 in
        p.plate_h_mm  = 482.6;   % 19 in

    case 'fine'
        % Color / mono board. 12 mm squares spans 76% x 66% of the color
        % camera FOV (174 x 145 mm at 0.071 mm/px, f=35 mm).
        % Cut from the 5 x 19 in offcut so the strip self-centres in the
        % 500 mm tank exactly as the 19 x 19 plate does.
        p.name        = 'fine';
        p.square_mm   = 12;
        p.n_cols      = 11;      % squares across  (odd)
        p.n_rows      = 8;       % squares down    (even)
        p.quiet_mm    = 10;      % narrower ring - strip is only 127 mm tall
        p.plate_w_mm  = 482.6;   % 19 in, along the tank
        p.plate_h_mm  = 127.0;   % 5 in, across the tank

    otherwise
        error('calib_target_params:badName', ...
            'Unknown board "%s" (expected ''coarse'' or ''fine'')', name);
end

%% Derived geometry
p.inner_cols = p.n_cols - 1;
p.inner_rows = p.n_rows - 1;
p.n_corners  = p.inner_cols * p.inner_rows;

p.board_w_mm = p.n_cols * p.square_mm;
p.board_h_mm = p.n_rows * p.square_mm;

% Board is centred on the plate; SVG origin is top-left, y down.
p.board_x0_mm = (p.plate_w_mm - p.board_w_mm) / 2;
p.board_y0_mm = (p.plate_h_mm - p.board_h_mm) / 2;

%% Inner corners in board-local coords (origin at board centre, mm)
[ix, iy] = meshgrid(1:p.inner_cols, 1:p.inner_rows);
p.corner_X = (ix(:) - (p.inner_cols + 1) / 2) * p.square_mm;
p.corner_Y = (iy(:) - (p.inner_rows + 1) / 2) * p.square_mm;

%% Validation
assert(mod(p.n_cols, 2) ~= mod(p.n_rows, 2), ...
    ['calib_target_params: %s board is %dx%d squares - must be odd x even ' ...
     'for unambiguous orientation.'], p.name, p.n_cols, p.n_rows);

margin_w = (p.plate_w_mm - p.board_w_mm) / 2 - p.quiet_mm;
margin_h = (p.plate_h_mm - p.board_h_mm) / 2 - p.quiet_mm;
assert(margin_w > 0 && margin_h > 0, ...
    ['calib_target_params: %s pattern + quiet ring (%.1f x %.1f mm) does ' ...
     'not fit the plate (%.1f x %.1f mm).'], p.name, ...
    p.board_w_mm + 2*p.quiet_mm, p.board_h_mm + 2*p.quiet_mm, ...
    p.plate_w_mm, p.plate_h_mm);
end
