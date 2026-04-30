function [Sx, Sy] = cisg_calculate_slopes(ref_image, obs_image, setup)
% CISG_CALCULATE_SLOPES Calculate surface slopes from CISG images
%
% Inputs:
%   ref_image - reference image (flat surface), RGB uint8 or uint16
%   obs_image - observed image (wavy surface), RGB uint8 or uint16
%   design - design parameters struct
%
% Outputs:
%   Sx - x-component of surface slope
%   Sy - y-component of surface slope

% Extract color channels and convert to double
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
mid_row = round(size(ref_image, 1) / 2);
R_row = double(ref_image(mid_row, :, 1));
cols = (1:length(R_row))';
p_x = polyfit(cols, R_row, 1);  % color units per pixel (x)

mid_col = round(size(ref_image, 2) / 2);
G_col = double(ref_image(:, mid_col, 2));
rows = (1:length(G_col))';
p_y = polyfit(rows, G_col, 1);  % color units per pixel (y)

% Convert using ruler calibration: (cm/pixel) / (color/pixel) = cm/color
scale_x = setup.cm_per_pixel_x / p_x(1);
scale_y = setup.cm_per_pixel_y / p_y(1);
% Physical displacement on pattern (cm)
delta_x = delta_R * scale_x;
delta_y = delta_G * scale_y;

% Convert displacement to slopes using refraction relationship
% delta = H * (n_water - n_air) / n_air * S
% Therefore: S = delta * n_air / (H * (n_water - n_air))

Sx = delta_x * setup.n_air / ...
    (setup.water_depth * (setup.n_water - setup.n_air));

Sy = delta_y * setup.n_air / ...
    (setup.water_depth * (setup.n_water - setup.n_air));

% fprintf('Slopes calculated from %d-bit images\n', log2(max_val+1));
% fprintf('  Sx range: [%.4f, %.4f]\n', min(Sx(:)), max(Sx(:)));
% fprintf('  Sy range: [%.4f, %.4f]\n', min(Sy(:)), max(Sy(:)));
% fprintf('  RMS slope magnitude: %.4f\n', std(sqrt(Sx(:).^2 + Sy(:).^2)));
end