%% calibrate_intrinsic.m
% Part 1 - Zhang's method on a checkerboard sequence recorded in air.
% Reads frames at manually selected indices from intrinsic.ats (or the
% synthetic intrinsic.mat for testing), runs detectCheckerboardPoints
% and estimateCameraParameters, prunes high-error poses, saves
% cameraParams to calibration/intrinsics.mat.

clear; clc;

%% ---- Config ----
cfg.pixel_pitch_mm            = 0.025;
cfg.nominal_focal_mm          = 13;
cfg.checker_inner_corners     = [7 5];   % [cols, rows]
cfg.checker_square_size_mm    = 30;
cfg.reproj_error_threshold_px = 1.0;

here = fileparts(mfilename('fullpath'));
cfg.intrinsic_path        = fullfile(here, '..', 'intrinsic.mat');   % swap to .ats in prod
cfg.intrinsic_frame_indices = 1:25;     % manual per-pose indices; edit for real data
cfg.out_dir               = here;
cfg.image_size            = [512 640];

addpath(here);   % so read_calib_frame is on path if run elsewhere

%% ---- Load + normalize frames ----
N = numel(cfg.intrinsic_frame_indices);
fprintf('Loading %d intrinsic frames from %s\n', N, cfg.intrinsic_path);

imgs_u8 = zeros(cfg.image_size(1), cfg.image_size(2), N, 'uint8');
for k = 1:N
    f = read_calib_frame(cfg.intrinsic_path, cfg.intrinsic_frame_indices(k));
    fmin = min(f(:)); fmax = max(f(:));
    if fmax > fmin
        imgs_u8(:,:,k) = uint8(255 * (f - fmin) / (fmax - fmin));
    else
        imgs_u8(:,:,k) = zeros(cfg.image_size, 'uint8');
    end
end

%% ---- Detect corners ----
fprintf('Detecting checkerboard corners...\n');
[imagePoints, boardSize, imagesUsed] = detectCheckerboardPoints( ...
    imgs_u8, 'PartialDetections', false);

expected_board = cfg.checker_inner_corners + 1;   % squares = inner + 1
if ~isequal(sort(boardSize-1), sort(cfg.checker_inner_corners))
    warning('Detected board %dx%d squares; expected inner corners %dx%d (i.e. %dx%d squares).', ...
        boardSize(1), boardSize(2), cfg.checker_inner_corners, expected_board);
end

drop_idx = find(~imagesUsed);
if ~isempty(drop_idx)
    fprintf('Detection failed on poses: %s\n', mat2str(cfg.intrinsic_frame_indices(drop_idx)));
end
fprintf('Detected on %d / %d poses\n', sum(imagesUsed), N);

%% ---- World points ----
worldPoints = generateCheckerboardPoints(boardSize, cfg.checker_square_size_mm);

%% ---- First fit ----
fprintf('Running estimateCameraParameters (first pass)...\n');
[cameraParams, ~, est_errs] = estimateCameraParameters( ...
    imagePoints, worldPoints, ...
    'ImageSize', cfg.image_size, ...
    'EstimateSkew', false, ...
    'NumRadialDistortionCoefficients', 2, ...
    'EstimateTangentialDistortion', false);

per_image_err = cameraParams.ReprojectionErrors;   % [Npts x 2 x Nimg]
mean_per_img  = squeeze(sqrt(mean(sum(per_image_err.^2, 2), 1)));

fprintf('\nPer-image reprojection error (px):\n');
used_pose_idx = cfg.intrinsic_frame_indices(imagesUsed);
for k = 1:numel(mean_per_img)
    fprintf('  pose %3d: %.3f px\n', used_pose_idx(k), mean_per_img(k));
end
fprintf('Mean error (first pass): %.3f px\n', cameraParams.MeanReprojectionError);

%% ---- Outlier pruning ----
keep = mean_per_img <= cfg.reproj_error_threshold_px;
if all(keep)
    fprintf('No poses exceed %.2f px threshold; keeping first-pass result.\n', ...
        cfg.reproj_error_threshold_px);
else
    fprintf('Dropping %d poses above %.2f px threshold; refitting...\n', ...
        sum(~keep), cfg.reproj_error_threshold_px);
    err_before = cameraParams.MeanReprojectionError;
    cameraParams = estimateCameraParameters( ...
        imagePoints(:,:,keep), worldPoints, ...
        'ImageSize', cfg.image_size, ...
        'EstimateSkew', false, ...
        'NumRadialDistortionCoefficients', 2, ...
        'EstimateTangentialDistortion', false);
    fprintf('Mean error: %.3f px -> %.3f px\n', err_before, cameraParams.MeanReprojectionError);
end

%% ---- Sanity prints ----
fl_px  = cameraParams.FocalLength;
fl_mm  = fl_px * cfg.pixel_pitch_mm;
pp_px  = cameraParams.PrincipalPoint;
pp_off = pp_px - (cfg.image_size([2 1]) + 1) / 2;
k_rad  = cameraParams.RadialDistortion;

fprintf('\n=== Intrinsic summary ===\n');
fprintf('Focal length: fx=%.3f px (%.3f mm)  fy=%.3f px (%.3f mm)  (nominal %d mm)\n', ...
    fl_px(1), fl_mm(1), fl_px(2), fl_mm(2), cfg.nominal_focal_mm);
fprintf('Principal point: (%.2f, %.2f) px; offset from center: (%+.2f, %+.2f) px\n', ...
    pp_px(1), pp_px(2), pp_off(1), pp_off(2));
fprintf('Radial distortion: [%.5f %.5f]\n', k_rad(1), k_rad(2));
fprintf('Mean reprojection error: %.3f px\n', cameraParams.MeanReprojectionError);

%% ---- Visualize ----
figure('Name','Reprojection errors');
showReprojectionErrors(cameraParams);
figure('Name','Extrinsics');
showExtrinsics(cameraParams, 'CameraCentric');

%% ---- Save ----
out_path = fullfile(cfg.out_dir, 'intrinsics.mat');
save(out_path, 'cameraParams', 'cfg', '-v7.3');
fprintf('\nSaved %s\n', out_path);
