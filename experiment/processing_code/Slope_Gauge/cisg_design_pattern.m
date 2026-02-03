function [design, pattern] = cisg_design_pattern(params)
% CISG_DESIGN_PATTERN Design CISG system and generate color pattern
%
% Inputs:
%   params - struct with fields:
%     .camera_height - cm above water surface
%     .water_depth - cm depth above pattern
%     .focal_length_mm - lens focal length
%     .pixel_width, pixel_height - camera resolution
%     .n_water, n_air - refractive indices
%
% Outputs:
%   design - struct with calculated design parameters
%   pattern - RGB image of CISG pattern

    % Camera sensor parameters (typical 1/2" sensor)
    sensor_width_mm = 8.8;
    sensor_height_mm = 6.6;
    
    % Calculate optical path and FOV
    total_height = params.camera_height + params.water_depth;
    
    % FOV in air
    FOV_width = 2 * total_height * tan(atan(sensor_width_mm/(2*params.focal_length_mm)));
    FOV_height = 2 * total_height * tan(atan(sensor_height_mm/(2*params.focal_length_mm)));
    
    % Refraction correction
    FOV_width_pattern = FOV_width * (1 + params.water_depth/params.camera_height * ...
        (params.n_water/params.n_air - 1));
    FOV_height_pattern = FOV_height * (1 + params.water_depth/params.camera_height * ...
        (params.n_water/params.n_air - 1));
    
    % Estimate max displacement (assume max slope ~0.3)
    max_slope = 0.3;
    max_displacement = total_height * tan(atan(max_slope)) * ...
        (params.n_water - params.n_air) / params.n_air;
    
    % Required pattern size
    pattern_width_cm = ceil((FOV_width_pattern + 2*max_displacement) / 5) * 5;
    pattern_height_cm = ceil((FOV_height_pattern + 2*max_displacement) / 5) * 5;
    
    % Spatial resolution
    spatial_res_x = FOV_width_pattern / params.pixel_width;
    spatial_res_y = FOV_height_pattern / params.pixel_height;
    
    % Print resolution
    desired_pattern_res = spatial_res_x / 4;
    dpi = max(300, ceil(2.54 / desired_pattern_res));
    
    % Store design parameters
    design = struct();
    design.camera_height = params.camera_height;
    design.water_depth = params.water_depth;
    design.FOV_width = FOV_width_pattern;
    design.FOV_height = FOV_height_pattern;
    design.pattern_width_cm = pattern_width_cm;
    design.pattern_height_cm = pattern_height_cm;
    design.spatial_res_x = spatial_res_x;
    design.spatial_res_y = spatial_res_y;
    design.n_water = params.n_water;
    design.n_air = params.n_air;
    design.pixel_width = params.pixel_width;
    design.pixel_height = params.pixel_height;
    design.dpi = dpi;
    
    % Print summary
    fprintf('Field of View: %.1f x %.1f cm\n', FOV_width_pattern, FOV_height_pattern);
    fprintf('Pattern size: %.1f x %.1f cm\n', pattern_width_cm, pattern_height_cm);
    fprintf('Spatial resolution: %.3f cm/pixel\n', spatial_res_x);
    fprintf('Recommended DPI: %d\n', dpi);
    
    % Generate pattern
    pattern_width_px = round(pattern_width_cm / 2.54 * dpi);
    pattern_height_px = round(pattern_height_cm / 2.54 * dpi);
    
    [X, Y] = meshgrid(1:pattern_width_px, 1:pattern_height_px);
    
    Red = uint8(255 * (X - 1) / (pattern_width_px - 1));
    Green = uint8(255 * (Y - 1) / (pattern_height_px - 1));
    Blue = uint8(128 * ones(pattern_height_px, pattern_width_px));
    
    pattern = cat(3, Red, Green, Blue);
    
    % Save pattern
    imwrite(pattern, 'CISG_pattern.png');
    fprintf('Pattern saved to CISG_pattern.png\n');
    
    % Save design parameters
    save('CISG_design.mat', 'design');
end
