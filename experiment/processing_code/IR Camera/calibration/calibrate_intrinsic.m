%% calibrate_intrinsic.m
% Part 1 - Zhang's method on a checkerboard sequence recorded in air.
% Reads frames at manually selected indices from intrinsic.ats (or the
% synthetic intrinsic.mat for testing), runs detectCheckerboardPoints
% and estimateCameraParameters, prunes high-error poses, saves
% cameraParams to calibration/intrinsics.mat.

clear; clc;

%% ---- Config ----
here = fileparts(mfilename('fullpath'));
addpath(here);   % so read_calib_frame / calib_target_params are on path

% Board geometry comes from calib_target_params.m, NOT hardcoded here -
% that is the single source of truth shared with the SVG cut file and
% corners.csv. If you change the physical target, this picks it up
% automatically; nothing below needs editing.
target = calib_target_params('coarse');

cfg.pixel_pitch_mm            = 0.025;
cfg.nominal_focal_mm          = 13;
cfg.checker_inner_corners     = [target.inner_cols, target.inner_rows];  % [cols, rows]
cfg.checker_square_size_mm    = target.square_mm;
cfg.reproj_error_threshold_px = 1.0;
cfg.num_radial_coeffs         = 3;   % visible fisheye-like bending in real
                                      % frames - the 2-coeff Brown-Conrady
                                      % truncation under-fits strong barrel
                                      % distortion at the periphery

cfg.intrinsic_path        = fullfile(here, '..', 'intrinsic.mat');   % swap to .ats in prod
cfg.intrinsic_frame_indices = 1:25;     % manual per-pose indices; edit for real data
cfg.out_dir               = here;
cfg.image_size            = [512 640];

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
    'NumRadialDistortionCoefficients', cfg.num_radial_coeffs, ...
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
        'NumRadialDistortionCoefficients', cfg.num_radial_coeffs, ...
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
fprintf('Radial distortion: [%s]\n', sprintf('%.5f ', k_rad));
fprintf('Mean reprojection error: %.3f px\n', cameraParams.MeanReprojectionError);

%% ---- Radial residual check ----
% A flat radial-error profile means the distortion model fits well
% everywhere. A profile that climbs with radius means the model is
% under-fitting the periphery - i.e. cfg.num_radial_coeffs is still too
% low for how strong the real (fisheye-like) distortion is, and points
% near the image edge should not yet be trusted.
imgPts_used  = imagePoints(:, :, keep);
reprojPts    = cameraParams.ReprojectedPoints;
errs_xy      = imgPts_used - reprojPts;
errs_px      = sqrt(sum(errs_xy.^2, 2));                    % [Npts x 1 x Nimg]
center       = (cfg.image_size([2 1]) + 1) / 2;
radius_px    = sqrt(sum((imgPts_used - reshape(center,1,2)).^2, 2));

radius_flat  = radius_px(:);
err_flat     = errs_px(:);
max_radius   = sqrt(sum((cfg.image_size([2 1]) / 2).^2));
edge_bins    = [0 0.5 0.75 1.0] * max_radius;
fprintf('\nRadial residual check (want roughly flat across bins):\n');
for b = 1:numel(edge_bins)-1
    in_bin = radius_flat >= edge_bins(b) & radius_flat < edge_bins(b+1);
    if any(in_bin)
        fprintf('  r in [%5.1f, %5.1f) px: mean err %.3f px  (n=%d)\n', ...
            edge_bins(b), edge_bins(b+1), mean(err_flat(in_bin)), sum(in_bin));
    end
end
outer = radius_flat >= edge_bins(1);
inner = radius_flat < edge_bins(2);
if mean(err_flat(~inner)) > 2 * mean(err_flat(inner))
    warning(['Edge reprojection error is >2x the center error even with ' ...
        '%d radial coefficients. The pinhole+Brown-Conrady model is ' ...
        'likely inadequate for this lens - consider estimateFisheyeParameters ' ...
        '/ the fisheyeParameters workflow instead.'], cfg.num_radial_coeffs);
end

%% ---- Visualize ----
figure('Name','Reprojection errors');
showReprojectionErrors(cameraParams);
figure('Name','Extrinsics');
showExtrinsics(cameraParams, 'CameraCentric');

%% ---- Save ----
out_path = fullfile(cfg.out_dir, 'intrinsics.mat');
save(out_path, 'cameraParams', 'cfg', '-v7.3');
fprintf('\nSaved %s\n', out_path);
