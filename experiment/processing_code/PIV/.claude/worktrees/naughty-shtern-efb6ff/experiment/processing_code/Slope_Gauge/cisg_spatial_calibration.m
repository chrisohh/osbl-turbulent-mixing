%% CISG_calibration.m - Pixel-to-cm calibration using ruler images
% clear all;

%% Setup
coreview_y = 6;   % ruler in y direction
coreview_x = 7;   % ruler in x direction

% Ruler units - change this if you switch rulers
ruler_unit = 'inches';  % 'inches', 'cm', 'mm'
switch ruler_unit
    case 'inches'
        unit_to_cm = 2.54;
    case 'cm'
        unit_to_cm = 1.0;
    case 'mm'
        unit_to_cm = 0.1;
end

%% Detect file naming format and load image (y-direction, CoreView 6)
ref_folderName_y = sprintf("\\\\airseaserver28\\D\\HLAB_2026\\SlopeGauge\\CoreView_%d\\Flare 12M125 CCL C1112A00\\",coreview_y);
files_y = dir(strcat(ref_folderName_y, sprintf('CoreView_%d_*.raw', coreview_y)));
first_file_y = files_y(1).name;
fprintf('Loading CoreView_%d (y-direction ruler)...\n', coreview_y);
[img_y, ~] = cisg_load_coreview(strcat(ref_folderName_y,first_file_y));

%% Pick two points on y-ruler
figure('Name', sprintf('CoreView_%d: Y-direction ruler - Click two points', coreview_y), ...
       'Position', [100 100 1200 800]);
imshow(img_y);
title(sprintf('CoreView_%d: Click two ruler points (y-direction)\nNote the ruler reading at each point', coreview_y));

fprintf('\nClick first point on y-ruler...\n');
[x1_y, y1_y] = ginput(1);
fprintf('  Point 1: pixel (%.1f, %.1f)\n', x1_y, y1_y);

fprintf('Click second point on y-ruler...\n');
[x2_y, y2_y] = ginput(1);
fprintf('  Point 2: pixel (%.1f, %.1f)\n', x2_y, y2_y);

% Enter ruler readings
ruler_val1_y = input(sprintf('Enter ruler reading at point 1 (%s): ', ruler_unit));
ruler_val2_y = input(sprintf('Enter ruler reading at point 2 (%s): ', ruler_unit));

% Calculate cm/pixel for y
dist_cm_y = abs(ruler_val2_y - ruler_val1_y) * unit_to_cm;
dist_pix_y = sqrt((x2_y - x1_y)^2 + (y2_y - y1_y)^2);
cm_per_pixel_y = dist_cm_y / dist_pix_y;

fprintf('\n--- Y-direction calibration ---\n');
fprintf('  Distance: %.2f %s = %.2f cm = %.1f pixels\n', ...
    abs(ruler_val2_y - ruler_val1_y), ruler_unit, dist_cm_y, dist_pix_y);
fprintf('  cm/pixel (y): %.4f\n\n', cm_per_pixel_y);

%% Load image (x-direction, CoreView 7)
ref_folderName_x = sprintf("\\\\airseaserver28\\D\\HLAB_2026\\SlopeGauge\\CoreView_%d\\Flare 12M125 CCL C1112A00\\", coreview_x);
files_x = dir(strcat(ref_folderName_x, sprintf('CoreView_%d_*.raw', coreview_x)));
first_file_x = files_x(1).name;

fprintf('Loading CoreView_%d (x-direction ruler)...\n', coreview_x);
[img_x, ~] = cisg_load_coreview(strcat(ref_folderName_x, first_file_x));

%% Pick two points on x-ruler
figure('Name', sprintf('CoreView_%d: X-direction ruler - Click two points', coreview_x), ...
       'Position', [100 100 1200 800]);
imshow(img_x);
title(sprintf('CoreView_%d: Click two ruler points (x-direction)\nNote the ruler reading at each point', coreview_x));

fprintf('Click first point on x-ruler...\n');
[x1_x, y1_x] = ginput(1);
fprintf('  Point 1: pixel (%.1f, %.1f)\n', x1_x, y1_x);

fprintf('Click second point on x-ruler...\n');
[x2_x, y2_x] = ginput(1);
fprintf('  Point 2: pixel (%.1f, %.1f)\n', x2_x, y2_x);

% Enter ruler readings
ruler_val1_x = input(sprintf('Enter ruler reading at point 1 (%s): ', ruler_unit));
ruler_val2_x = input(sprintf('Enter ruler reading at point 2 (%s): ', ruler_unit));

% Calculate cm/pixel for x
dist_cm_x = abs(ruler_val2_x - ruler_val1_x) * unit_to_cm;
dist_pix_x = sqrt((x2_x - x1_x)^2 + (y2_x - y1_x)^2);
cm_per_pixel_x = dist_cm_x / dist_pix_x;

fprintf('\n--- X-direction calibration ---\n');
fprintf('  Distance: %.2f %s = %.2f cm = %.1f pixels\n', ...
    abs(ruler_val2_x - ruler_val1_x), ruler_unit, dist_cm_x, dist_pix_x);
fprintf('  cm/pixel (x): %.4f\n\n', cm_per_pixel_x);

%% Summary and save
fprintf('===== CALIBRATION SUMMARY =====\n');
fprintf('  cm/pixel (x): %.4f\n', cm_per_pixel_x);
fprintf('  cm/pixel (y): %.4f\n', cm_per_pixel_y);
fprintf('  Ratio (x/y):  %.4f (should be ~1 if pixels are square)\n\n', cm_per_pixel_x/cm_per_pixel_y);

cal = struct();
cal.cm_per_pixel_x = cm_per_pixel_x;
cal.cm_per_pixel_y = cm_per_pixel_y;
cal.coreview_x = coreview_x;
cal.coreview_y = coreview_y;
cal.ruler_unit = ruler_unit;
cal.unit_to_cm = unit_to_cm;
% Store raw picking info for reference
cal.points_x = struct('pix', [x1_x y1_x; x2_x y2_x], 'ruler_vals', [ruler_val1_x ruler_val2_x]);
cal.points_y = struct('pix', [x1_y y1_y; x2_y y2_y], 'ruler_vals', [ruler_val1_y ruler_val2_y]);

save(sprintf('CISG_Calibration_CoreView_%d.mat', coreview_y), 'cal', '-v7.3');
fprintf('Saved: CISG_calibration.mat\n');

%% Visualize picked points
figure('Name', 'Calibration Verification', 'Position', [100 100 1600 700]);

subplot(1,2,1);
imshow(img_y);
hold on;
plot([x1_y x2_y], [y1_y y2_y], 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
title(sprintf('CoreView_%d (y): %.4f cm/pixel', coreview_y, cm_per_pixel_y));

subplot(1,2,2);
imshow(img_x);
hold on;
plot([x1_x x2_x], [y1_x y2_x], 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
title(sprintf('CoreView_%d (x): %.4f cm/pixel', coreview_x, cm_per_pixel_x));

fprintf('\n✓ Calibration complete!\n');