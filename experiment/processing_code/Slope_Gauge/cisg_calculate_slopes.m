function [Sx, Sy, info] = cisg_calculate_slopes(ref_image, obs_image, setup)
% CISG_CALCULATE_SLOPES Calculate surface slopes from CISG images
%
% Inputs:
%   ref_image - reference image (flat surface), RGB uint8 or uint16
%   obs_image - observed image (wavy surface), RGB uint8 or uint16
%   design - design parameters struct
%
% Outputs (info: colour-inversion info of the 'local' map, e.g. info.dye_c with setup.dye_model = 'joint'; [] otherwise):
%   Sx - d eta / dx, x = image right (columns, downwind)
%   Sy - d eta / dy, y = image UP (the tank frame of the plate calibration).  The colour shift is measured along
%        rows (down), so Sy = -(row shift) x scale.  Before 2026-10-25 Sy had the opposite sign (d eta / d row):
%        the curl of (Sx, Sy) was ~1 x its divergence and ~0.3 with Sy flipped (cisg_pipeline 5j); caches made
%        before then are not reused (folder tag _yup, flat-bias key sy_up).

% Extract color channels and convert to double
info = [];
R_ref = double(ref_image(:,:,1));
G_ref = double(ref_image(:,:,2));

R_obs = double(obs_image(:,:,1));
G_obs = double(obs_image(:,:,2));

% Determine max value based on image data type
if isa(ref_image, 'uint16')
    max_val = 65535;
else
    max_val = 255;
end

% Calculate color displacement (in units of max_val)
delta_R = R_obs - R_ref;
delta_G = G_obs - G_ref;

% Calibration: color value to physical displacement
% Red channel spans 0-max_val over pattern width
% Green channel spans 0-max_val over pattern height
% scale_x = design.pattern_width_cm / max_val;   % cm per color unit
% scale_y = design.pattern_height_cm / max_val;  % cm per color unit

% For the Red channel (x-direction):
% Grab a row near the center of the image (avoid edges)
% Calibration: color value to physical displacement
% Fit color gradient from reference image
use_new = ~isfield(setup, 'use_apparent_depth') || setup.use_apparent_depth;
if use_new && ~isfield(setup, 'L_eff_cm')
    error('cisg_calculate_slopes:L_eff', ...
        'setup.use_apparent_depth needs setup.L_eff_cm (cisg_apparent_depth); set it false for the legacy path.');
end

if use_new && isfield(setup, 'colour_map') && strcmp(setup.colour_map, 'local')
    % Local colour -> position inversion (Jacobian of the reference frame); NaN where untrustworthy.
    % dpx: shift along columns, dpy: along rows (same sign as the linear path), in pixels
    if isfield(setup, 'init_Sx') && ~isempty(setup.init_Sx)        % warm start: previous frame's slopes (ROI-sized) -> shift (px)
        setup.init_dx = double(setup.init_Sx) * ((setup.n_water - setup.n_air) * setup.L_eff_cm) / (setup.pattern_scale * setup.cm_per_pixel_x);
        setup.init_dy = -double(setup.init_Sy) * ((setup.n_water - setup.n_air) * setup.L_eff_cm) / (setup.pattern_scale * setup.cm_per_pixel_y);
    end
    [dpx, dpy, info] = cisg_colour_displacement(ref_image, obs_image, setup);
    % the ray bends uphill, so the pattern seen at x moves by +k grad(eta): columns run with +x, rows against +y
    Sx =  dpx * (setup.pattern_scale * setup.cm_per_pixel_x) / ((setup.n_water - setup.n_air) * setup.L_eff_cm);
    Sy = -dpy * (setup.pattern_scale * setup.cm_per_pixel_y) / ((setup.n_water - setup.n_air) * setup.L_eff_cm);
    k_sp = 4;  if isfield(setup, 'speckle_k'), k_sp = setup.speckle_k; end        % isolated-speckle filter (0 = off)
    info.n_speckle = 0;
    if k_sp > 0
        rr = [1 size(Sx, 1)];  cc = [1 size(Sx, 2)];
        if isfield(setup, 'roi'), rr = [setup.roi.row_min setup.roi.row_max];  cc = [setup.roi.col_min setup.roi.col_max]; end
        [a, n1] = cisg_despeckle(Sx(rr(1):rr(2), cc(1):cc(2)), k_sp);  Sx(rr(1):rr(2), cc(1):cc(2)) = a;
        [a, n2] = cisg_despeckle(Sy(rr(1):rr(2), cc(1):cc(2)), k_sp);  Sy(rr(1):rr(2), cc(1):cc(2)) = a;
        info.n_speckle = [n1 n2];
    end
    return;
end

mid_row = round(size(ref_image, 1) / 2);
R_row = double(ref_image(mid_row, :, 1));
cols = (1:length(R_row))';
fit_cols = true(size(cols));
if use_new && isfield(setup, 'target_cols')
    % fit only inside the target (inset 5%): the surround is not part of the ramp
    tc = setup.target_cols;  inset = 0.05 * diff(tc);
    fit_cols = cols >= tc(1) + inset & cols <= tc(2) - inset;
end
p_x = polyfit(cols(fit_cols), R_row(fit_cols)', 1);  % color units per pixel (x)

mid_col = round(size(ref_image, 2) / 2);
G_col = double(ref_image(:, mid_col, 2));
rows = (1:length(G_col))';
p_y = polyfit(rows, G_col, 1);  % color units per pixel (y)

if use_new
    % pixel size ON THE PATTERN = surface cm/px * r (the target is deeper than the surface)
    r = setup.pattern_scale;
else
    r = 1;
end
% Convert using ruler calibration: (cm/pixel) / (color/pixel) = cm/color
scale_x = r * setup.cm_per_pixel_x / p_x(1);
scale_y = r * setup.cm_per_pixel_y / p_y(1);
% Physical displacement on pattern (cm); the green ramp is fitted along rows (down), so negate it for y up
delta_x = delta_R * scale_x;
delta_y = -delta_G * scale_y;

if use_new
    % A ray is deviated by S (n_w - n_a)/n_w in the water, so the pattern is displaced by
    %   delta = S (n_w - n_a) L_eff,   L_eff = h/n_w + t/n_g + g  (air-equivalent depth)
    % hence S = delta / ((n_w - n_a) L_eff)   [cisg_apparent_depth.m]
    Sx = delta_x / ((setup.n_water - setup.n_air) * setup.L_eff_cm);
    Sy = delta_y / ((setup.n_water - setup.n_air) * setup.L_eff_cm);
else
    % Legacy: delta = H * (n_water - n_air) / n_air * S
    % Therefore: S = delta * n_air / (H * (n_water - n_air))
    Sx = delta_x * setup.n_air / ...
        (setup.water_depth * (setup.n_water - setup.n_air));

    Sy = delta_y * setup.n_air / ...
        (setup.water_depth * (setup.n_water - setup.n_air));
end

% fprintf('Slopes calculated from %d-bit images\n', log2(max_val+1));
% fprintf('  Sx range: [%.4f, %.4f]\n', min(Sx(:)), max(Sx(:)));
% fprintf('  Sy range: [%.4f, %.4f]\n', min(Sy(:)), max(Sy(:)));
% fprintf('  RMS slope magnitude: %.4f\n', std(sqrt(Sx(:).^2 + Sy(:).^2)));
end