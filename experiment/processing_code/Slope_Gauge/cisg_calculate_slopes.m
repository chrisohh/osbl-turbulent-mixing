function [Sx, Sy] = cisg_calculate_slopes(ref_image, obs_image, design)
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
scale_x = design.pattern_width_cm / max_val;   % cm per color unit
scale_y = design.pattern_height_cm / max_val;  % cm per color unit

% Physical displacement on pattern (cm)
delta_x = delta_R * scale_x;
delta_y = delta_G * scale_y;

% Convert displacement to slopes using refraction relationship
% delta = H * (n_water - n_air) / n_air * S
% Therefore: S = delta * n_air / (H * (n_water - n_air))

Sx = delta_x * design.n_air / ...
    (design.water_depth * (design.n_water - design.n_air));

Sy = delta_y * design.n_air / ...
    (design.water_depth * (design.n_water - design.n_air));

fprintf('Slopes calculated from %d-bit images\n', log2(max_val+1));
fprintf('  Sx range: [%.4f, %.4f]\n', min(Sx(:)), max(Sx(:)));
fprintf('  Sy range: [%.4f, %.4f]\n', min(Sy(:)), max(Sy(:)));
fprintf('  RMS slope magnitude: %.4f\n', std(sqrt(Sx(:).^2 + Sy(:).^2)));
end