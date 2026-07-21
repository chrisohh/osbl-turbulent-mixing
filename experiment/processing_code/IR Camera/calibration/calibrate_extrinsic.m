%% calibrate_extrinsic.m
% Part 2 - planar homography from undistorted image corners to tank-frame
% physical coordinates. Requires intrinsics.mat from calibrate_intrinsic.m.

clear; clc;

%% ---- Config ----
here = fileparts(mfilename('fullpath'));
cfg.intrinsics_path     = fullfile(here, 'intrinsics.mat');
cfg.extrinsic_path      = fullfile(here, '..', 'extrinsic.mat');   % swap to .ats in prod
cfg.extrinsic_frame_idx = 1;
cfg.corners_csv         = fullfile(here, '..', 'corners.csv');
cfg.out_dir             = here;
cfg.image_size          = [512 640];
cfg.grid_spacing_mm     = 50;   % overlay grid spacing for verification

addpath(here);

%% ---- Load intrinsics ----
S = load(cfg.intrinsics_path, 'cameraParams');
cameraParams = S.cameraParams;

%% ---- Load + normalize extrinsic frame ----
frame = read_calib_frame(cfg.extrinsic_path, cfg.extrinsic_frame_idx);
fmin = min(frame(:)); fmax = max(frame(:));
frame_u8 = uint8(255 * (frame - fmin) / max(fmax - fmin, eps));

%% ---- Detect corners ----
[imagePoints, boardSize] = detectCheckerboardPoints(frame_u8, 'PartialDetections', false);
expected = prod(boardSize - 1);
if size(imagePoints,1) ~= expected || any(isnan(imagePoints(:)))
    error('calibrate_extrinsic:detection', ...
        'Expected %d corners; detected %d. Inspect extrinsic frame.', ...
        expected, size(imagePoints,1));
end
fprintf('Detected %d corners (board %dx%d squares)\n', ...
    size(imagePoints,1), boardSize(1), boardSize(2));

%% ---- Undistort (mandatory before fitgeotrans) ----
undistortedPoints = undistortPoints(imagePoints, cameraParams);

%% ---- Read tank-frame physical corners ----
T = readtable(cfg.corners_csv);
if ~all(ismember({'X_mm','Y_mm'}, T.Properties.VariableNames))
    error('corners.csv must have headers X_mm,Y_mm');
end
physical_xy = [T.X_mm, T.Y_mm];
if size(physical_xy,1) ~= size(undistortedPoints,1)
    error('corners.csv has %d rows but %d corners were detected', ...
        size(physical_xy,1), size(undistortedPoints,1));
end

%% ---- Fit homography ----
H = fitgeotrans(undistortedPoints, physical_xy, 'projective');

%% ---- Validate ----
xy_pred = transformPointsForward(H, undistortedPoints);
residuals = xy_pred - physical_xy;
rms_mm = sqrt(mean(sum(residuals.^2, 2)));

fprintf('\nHomography fit:\n');
fprintf('  RMS residual: %.3f mm (target < 1 mm)\n', rms_mm);
if rms_mm > 1
    warning('RMS residual %.3f mm exceeds 1 mm target.', rms_mm);
end

% Pixel-to-mm scale at image center, for a quick sanity check
center_px = (cfg.image_size([2 1]) + 1) / 2;
center_px_u = undistortPoints(center_px, cameraParams);
nearby_u    = undistortPoints(center_px + [1 0], cameraParams);
near_mm     = transformPointsForward(H, [center_px_u; nearby_u]);
scale_mm_per_px_center = norm(near_mm(2,:) - near_mm(1,:));
fprintf('  Pixel-to-mm scale at image center: %.3f mm/px\n', scale_mm_per_px_center);

%% ---- Visualize: regular tank-frame grid back-projected to image ----
x_range = [min(physical_xy(:,1)) max(physical_xy(:,1))];
y_range = [min(physical_xy(:,2)) max(physical_xy(:,2))];
% expand by one grid spacing on each side
x_range = x_range + cfg.grid_spacing_mm * [-1 1];
y_range = y_range + cfg.grid_spacing_mm * [-1 1];

gx = floor(x_range(1)/cfg.grid_spacing_mm)*cfg.grid_spacing_mm : ...
     cfg.grid_spacing_mm : ...
     ceil(x_range(2)/cfg.grid_spacing_mm)*cfg.grid_spacing_mm;
gy = floor(y_range(1)/cfg.grid_spacing_mm)*cfg.grid_spacing_mm : ...
     cfg.grid_spacing_mm : ...
     ceil(y_range(2)/cfg.grid_spacing_mm)*cfg.grid_spacing_mm;

figure('Name','Extrinsic verification');
imagesc(frame); colormap(hot); axis image; hold on;
title(sprintf('Tank-frame %d mm grid back-projected to image (rms=%.2f mm)', ...
    cfg.grid_spacing_mm, rms_mm));

% Horizontal lines (constant y, varying x)
for yy = gy
    pts_phys = [gx(:), repmat(yy, numel(gx), 1)];
    pts_undist = transformPointsInverse(H, pts_phys);
    % re-apply distortion to land on raw image coords
    pts_img = redistortPoints(pts_undist, cameraParams);
    plot(pts_img(:,1), pts_img(:,2), 'c-', 'LineWidth', 0.5);
end
% Vertical lines (constant x, varying y)
for xx = gx
    pts_phys = [repmat(xx, numel(gy), 1), gy(:)];
    pts_undist = transformPointsInverse(H, pts_phys);
    pts_img = redistortPoints(pts_undist, cameraParams);
    plot(pts_img(:,1), pts_img(:,2), 'c-', 'LineWidth', 0.5);
end
plot(imagePoints(:,1), imagePoints(:,2), 'g+', 'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('Pixel u'); ylabel('Pixel v');

%% ---- Save ----
out_path = fullfile(cfg.out_dir, 'homography.mat');
save(out_path, 'H', 'rms_mm', 'undistortedPoints', 'physical_xy', ...
    'imagePoints', 'cfg', '-v7.3');
fprintf('Saved %s\n', out_path);

%% ---- helper ----
function pts_dist = redistortPoints(pts_undist, cameraParams)
% Forward-distort already-undistorted pixel points using the camera's
% radial+tangential model. MATLAB exposes undistortPoints but not the
% forward direction; we apply the standard pinhole model directly.
    K  = cameraParams.IntrinsicMatrix';
    fx = K(1,1); fy = K(2,2); cx = K(1,3); cy = K(2,3);
    kr = cameraParams.RadialDistortion;
    if numel(kr) < 3, kr(3) = 0; end
    if isprop(cameraParams,'TangentialDistortion')
        p = cameraParams.TangentialDistortion;
    else
        p = [0 0];
    end
    x = (pts_undist(:,1) - cx) / fx;
    y = (pts_undist(:,2) - cy) / fy;
    r2 = x.^2 + y.^2;
    radial = 1 + kr(1)*r2 + kr(2)*r2.^2 + kr(3)*r2.^3;
    x_d = x .* radial + 2*p(1).*x.*y + p(2)*(r2 + 2*x.^2);
    y_d = y .* radial + p(1)*(r2 + 2*y.^2) + 2*p(2).*x.*y;
    pts_dist = [x_d*fx + cx, y_d*fy + cy];
end
