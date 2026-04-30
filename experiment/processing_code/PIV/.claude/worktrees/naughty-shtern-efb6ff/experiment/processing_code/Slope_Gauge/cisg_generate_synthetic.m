function [ref_image, obs_image, true_slopes] = cisg_generate_synthetic(design, pattern, wave_params)
% CISG_GENERATE_SYNTHETIC Generate synthetic wave data for testing
%
% Inputs:
%   design - design parameters from cisg_design_pattern
%   pattern - RGB pattern image
%   wave_params - struct with:
%     .amplitude - wave amplitude (cm)
%     .wavelength_x, wavelength_y - wavelengths (cm)
%     .angle - wave angle (degrees)
%
% Outputs:
%   ref_image - reference image (flat surface)
%   obs_image - observed image (wavy surface)
%   true_slopes - struct with .Sx, .Sy, .eta, .x, .y

    % Create wave field
    x = linspace(0, design.FOV_width, design.pixel_width);
    y = linspace(0, design.FOV_height, design.pixel_height);
    [X, Y] = meshgrid(x, y);
    
    % Wave numbers
    kx = 2*pi / wave_params.wavelength_x;
    ky = 2*pi / wave_params.wavelength_y;
    
    % Multi-component wave field
    eta = wave_params.amplitude * (...
        sin(kx * X) + ...
        0.7 * sin(ky * Y) + ...
        0.5 * sin(kx * cosd(wave_params.angle) * X + ...
                   ky * sind(wave_params.angle) * Y) + ...
        0.3 * sin(2*kx * X + 1.5*ky * Y));
    
    % Calculate true slopes
    [Sx_true, Sy_true] = gradient(eta);
    Sx_true = Sx_true / mean(diff(x));
    Sy_true = Sy_true / mean(diff(y));
    
    fprintf('Wave RMS elevation: %.3f cm\n', std(eta(:)));
    fprintf('Wave RMS slope: %.4f\n', std(sqrt(Sx_true(:).^2 + Sy_true(:).^2)));
    
    % Store ground truth
    true_slopes = struct();
    true_slopes.Sx = Sx_true;
    true_slopes.Sy = Sy_true;
    true_slopes.eta = eta;
    true_slopes.x = x;
    true_slopes.y = y;
    
    % Generate reference image (flat surface)
    ref_image = sample_pattern(design, pattern, X, Y, 0, 0);
    
    % Calculate refraction displacement
    delta_x = design.water_depth * (design.n_water - design.n_air) / ...
              design.n_air * Sx_true;
    delta_y = design.water_depth * (design.n_water - design.n_air) / ...
              design.n_air * Sy_true;
    
    % Generate observed image (wavy surface)
    obs_image = sample_pattern(design, pattern, X + delta_x, Y + delta_y, 0, 0);
    
    fprintf('Synthetic images generated\n');
end

function img = sample_pattern(design, pattern, X_pattern, Y_pattern, offset_x, offset_y)
    % Sample pattern at given positions using bilinear interpolation
    
    [pattern_h, pattern_w, ~] = size(pattern);
    img = zeros(size(X_pattern, 1), size(X_pattern, 2), 3, 'uint8');
    
    % Convert physical coordinates to pixel coordinates
    % Pattern is centered on FOV
    X_px = (X_pattern + offset_x + design.pattern_width_cm/2) / ...
           design.pattern_width_cm * pattern_w;
    Y_px = (Y_pattern + offset_y + design.pattern_height_cm/2) / ...
           design.pattern_height_cm * pattern_h;
    
    % Bilinear interpolation
    for i = 1:size(img, 1)
        for j = 1:size(img, 2)
            px = X_px(i,j);
            py = Y_px(i,j);
            
            % Check bounds
            if px >= 1 && px <= pattern_w && py >= 1 && py <= pattern_h
                px1 = floor(px); px2 = min(ceil(px), pattern_w);
                py1 = floor(py); py2 = min(ceil(py), pattern_h);
                
                wx = px - px1;
                wy = py - py1;
                
                for c = 1:3
                    img(i,j,c) = uint8(...
                        (1-wx)*(1-wy)*double(pattern(py1,px1,c)) + ...
                        wx*(1-wy)*double(pattern(py1,px2,c)) + ...
                        (1-wx)*wy*double(pattern(py2,px1,c)) + ...
                        wx*wy*double(pattern(py2,px2,c)));
                end
            end
        end
    end
end
