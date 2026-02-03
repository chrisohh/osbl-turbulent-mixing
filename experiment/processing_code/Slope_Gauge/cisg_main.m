%% CISG Main Script - Complete Workflow for CoreView Data
clear all;

%% Setup parameters
coreview_num = 22;
coreview_ref_num = 19;%coreview_num
ref_fileName= sprintf("\\\\airseaserver28\\D\\HLAB_2026\\SlopeGauge\\CoreView_%d\\Flare 12M125 CCL C1112A00\\CoreView_%d_Flare 12M125 CCL C1112A00_01.raw", coreview_ref_num,coreview_ref_num);

obs_folderName = sprintf("\\\\airseaserver28\\D\\HLAB_2026\\SlopeGauge\\CoreView_%d\\Flare 12M125 CCL C1112A00\\", coreview_num);
fs = 50;  % Hz
dt = 1/fs;  % 0.02 seconds between frames

% Physical setup (ADJUST THESE TO YOUR ACTUAL MEASUREMENTS)
setup = struct();
setup.camera_height = 41.7;            % cm above water surface (MEASURE THIS)
setup.water_depth = 49;              % cm from water surface to pattern (MEASURE THIS)
% setup.pattern_width_cm = 20*2.54;    % 50.8 cm (your printed pattern)
% setup.pattern_height_cm = 20*2.54;   % 50.8 cm
setup.n_water = 1.333;
setup.n_air = 1.000;
setup.frame_rate = fs;
setup.dt = dt;

fprintf('===== CISG PROCESSING SETUP =====\n');
fprintf('Data folder: %s\n', ref_folderName);
fprintf('Frame rate: %d Hz (dt = %.4f s)\n', fs, dt);
fprintf('Camera height: %.1f cm\n', setup.camera_height);
fprintf('Water depth: %.1f cm\n', setup.water_depth);
% fprintf('Pattern size: %.1f x %.1f cm\n\n', setup.pattern_width_cm, setup.pattern_height_cm);

%% Load calibration
% mat produced from cisg_spatial_calibration.m
cal = load('CISG_Calibration_CoreView_6.mat');
cal = cal.cal;
setup.cm_per_pixel_x = cal.cm_per_pixel_x;
setup.cm_per_pixel_y = cal.cm_per_pixel_y;
%% Load and display reference image
fprintf('Loading reference image: %s\n', ref_fileName);
[ref_image, meta_ref] = cisg_load_coreview(ref_fileName);
[img_height, img_width] = size(ref_image, [1 2]);

% Display reference
figure('Name', 'Reference Image', 'Position', [100 100 1200 800]);
subplot(2,2,1);
imshow(ref_image);
title('Full RGB Reference');

subplot(2,2,2);
imagesc(ref_image(:,:,1));
colorbar;
title('Red Channel (increases left→right)');
axis equal tight;
colormap(gca, 'gray');

subplot(2,2,3);
imagesc(ref_image(:,:,2));
colorbar;
title('Green Channel (increases bottom→top)');
axis equal tight;
colormap(gca, 'gray');

subplot(2,2,4);
imagesc(ref_image(:,:,3));
colorbar;
title('Blue Channel (constant ~128)');
axis equal tight;
colormap(gca, 'gray');

fprintf('Reference image loaded: %d x %d\n\n', size(ref_image,1), size(ref_image,2));

%% find where to crop from the reference image
% Define valid region - click 4 corners on reference
figure('Name', 'Select valid region - Click 4 corners', 'Position', [100 100 1200 800]);
imshow(ref_image);
title('Click 4 corners of valid region (any order)');

fprintf('Click 4 corners...\n');
[cx, cy] = ginput(4);

% Sort points into left/right and top/bottom
% Left side: 2 points with smallest x, right side: 2 with largest x
[cx_sorted, idx] = sort(cx);
left_pts = idx(1:2);   % two leftmost
right_pts = idx(3:4);  % two rightmost

col_min = ceil(max(cx(left_pts)));    % rightmost of the left points
col_max = floor(min(cx(right_pts)));  % leftmost of the right points

[cy_sorted, idx] = sort(cy);
top_pts = idx(1:2);    % two smallest y (top of image)
bottom_pts = idx(3:4); % two largest y (bottom of image)

row_min = ceil(max(cy(top_pts)));     % lowest of the top points
row_max = floor(min(cy(bottom_pts))); % highest of the bottom points

setup.roi = struct('row_min', row_min, 'row_max', row_max, ...
    'col_min', col_min, 'col_max', col_max);

% Create mask
setup.valid_mask = false(img_height, img_width);
setup.valid_mask(row_min:row_max, col_min:col_max) = true;

fprintf('  ROI: rows [%d, %d], cols [%d, %d]\n', row_min, row_max, col_min, col_max);

%% Process frames
totalFrames = 2291;
start_frame = 1421;
end_frame = totalFrames;
frames_to_process = start_frame:10:end_frame;
nFrames = length(frames_to_process);

time = (frames_to_process - 1) * dt;

fprintf('===== PROCESSING FRAMES %d to %d (%d frames) =====\n', start_frame, end_frame, nFrames);

% Write directly to disk
% mat_fileName = sprintf('CISG_slopes_CoreView_%d.mat', coreview_num);
% m = matfile(mat_fileName, 'Writable', true);
% initialize zeros with 2 pages to force it to recognize 3D matrix
% m.Sx = zeros(img_height, img_width, 2, 'single');
% m.Sy = zeros(img_height, img_width, 2, 'single');

tic;
for count = 1:nFrames
    n = frames_to_process(count);
    obs_fileName = sprintf("CoreView_%d_Flare 12M125 CCL C1112A00_%04d.raw", coreview_num, n);
    [obs_image, ~] = cisg_load_coreview(strcat(sprintf("\\\\airseaserver28\\D\\HLAB_2026\\SlopeGauge\\CoreView_%d\\CoreView_%d\\Flare 12M125 CCL C1112A00\\", coreview_num,coreview_num), obs_fileName));
    
    [sx, sy] = cisg_calculate_slopes(ref_image, obs_image, setup);

    m.Sx(:,:,count+142) = single(sx);
    m.Sy(:,:,count+142) = single(sy);
   
    if mod(count, 10) == 0
        elapsed = toc;
        fps_processing = count / elapsed;
        est_remaining = (nFrames - count) / fps_processing;
        fprintf('  Frame %d/%d (%.1f%%) - %.1f fps - ETA: %.1f sec\n', ...
            count, nFrames, 100*count/nFrames, fps_processing, est_remaining);
    end
end
total_time = toc;

m.time = time;
m.setup = setup;

fprintf('\n✓ Processing complete!\n');
fprintf('  Total time: %.1f seconds\n', total_time);
fprintf('  Average speed: %.1f fps\n', nFrames/total_time);
fprintf('  Saved: %s\n\n', mat_fileName);

% Process each frame
% tic;
% count=0;
% for n = 1:totalFrames
%     count=1+count;
%     % Load observed image
%     obs_fileName = sprintf("CoreView_22_Flare 12M125 CCL C1112A00_%04d.raw", n);
%     [obs_image, meta_obs] = cisg_load_coreview(strcat(obs_folderName, obs_fileName));
% 
%     % Calculate slopes
%     [Sx(:,:,count), Sy(:,:,count)] = cisg_calculate_slopes(ref_image, obs_image, setup);
% 
%     % Progress update every 50 frames
%     if mod(n, 50) == 0
%         elapsed = toc;
%         fps_processing = n / elapsed;
%         est_remaining = (totalFrames - n) / fps_processing;
%         fprintf('  Frame %d/%d (%.1f%%) - %.1f fps - ETA: %.1f sec\n', ...
%             n, totalFrames, 100*n/totalFrames, fps_processing, est_remaining);
%     end
% end
% total_time = toc;

% fprintf('\n✓ Processing complete!\n');
% fprintf('  Total time: %.1f seconds\n', total_time);
% fprintf('  Average speed: %.1f fps\n', totalFrames/total_time);

%% Save results
% fprintf('\n===== SAVING RESULTS =====\n');
% save('CISG_slopes_frame1to1415_increment10.mat', 'Sx', 'Sy', 'time', 'setup', '-v7.3');
% fprintf('Saved: CISG_slopes_all_frames.mat\n');
% 
% % Save metadata
% fid = fopen('CISG_processing_info.txt', 'w');
% fprintf(fid, 'CISG Processing Summary\n');
% fprintf(fid, '======================\n\n');
% fprintf(fid, 'Data folder: %s\n', ref_folderName);
% fprintf(fid, 'Reference file: %s\n', ref_fileName);
% fprintf(fid, 'Total frames: %d\n', totalFrames);
% fprintf(fid, 'Frame rate: %d Hz\n', fs);
% fprintf(fid, 'Duration: %.2f seconds\n', totalFrames * dt);
% fprintf(fid, 'Image size: %d x %d\n', height, width);
% fprintf(fid, 'Camera height: %.1f cm\n', setup.camera_height);
% fprintf(fid, 'Water depth: %.1f cm\n', setup.water_depth);
% % fprintf(fid, 'Pattern size: %.1f x %.1f cm\n', setup.pattern_width_cm, setup.pattern_height_cm);
% fprintf(fid, 'Processing time: %.1f seconds\n', total_time);
% fclose(fid);
% fprintf('Saved: CISG_processing_info.txt\n');

%% Quick visualization of results
fprintf('\n===== VISUALIZING SAMPLE FRAMES =====\n');

figure('Name', 'Spatial Maps', 'Position', [100 100 1600 1000]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

frames_to_show = [228];%[1, round(size(m.Sx,3)/2), size(m.Sx,3)];

% Physical axes centered on ROI center
nx = setup.roi.col_max - setup.roi.col_min + 1;
ny = setup.roi.row_max - setup.roi.row_min + 1;
x_cm = ((1:nx) - nx/2) * setup.cm_per_pixel_x;
y_cm = ((1:ny) - ny/2) * setup.cm_per_pixel_y;

for i = 1:length(frames_to_show)
    Sx_frame = double(m.Sx(setup.roi.row_min:setup.roi.row_max, ...
                           setup.roi.col_min:setup.roi.col_max, frames_to_show(i)));
    
    nexttile;
    imagesc(x_cm, y_cm, Sx_frame);
    axis equal tight; colorbar;
    title(sprintf('S_x: Frame %d (t=%.2fs)', frames_to_show(i)*10, time(frames_to_show(i))));
    colormap(gca, 'turbo');
    caxis([-0.2,0.2]);
end

for i = 1:length(frames_to_show)
    Sy_frame = double(m.Sy(setup.roi.row_min:setup.roi.row_max, ...
                           setup.roi.col_min:setup.roi.col_max, frames_to_show(i)));
    
    nexttile;
    imagesc(x_cm, y_cm, Sy_frame);
    axis equal tight; colorbar;
    title(sprintf('S_y: Frame %d (t=%.2fs)', frames_to_show(i)*10, time(frames_to_show(i))));
    colormap(gca, 'turbo');
    caxis([-0.2,0.2]);
end


for i = 1:length(frames_to_show)
    frame_num = frames_to_show(i);
    Smag = sqrt(double(m.Sx(:,:,frames_to_show(i))).^2 + double(m.Sy(:,:,frames_to_show(i))).^2);
    Smag = Smag(setup.roi.row_min:setup.roi.row_max, ...
                           setup.roi.col_min:setup.roi.col_max);
    
    nexttile;
    imagesc(x_cm, y_cm, Smag);
    axis equal tight; colorbar;
    title(sprintf('|S|: Frame %d (t=%.2fs)', frames_to_show(i)*10, time(frames_to_show(i))));
    colormap(gca, 'hot');
    caxis([0,0.3]);
end

%% Make a video
opts = struct();
opts.caxis_limits = [-0.2, 0.2];      % Sx/Sy limits
opts.mag_caxis_limits = [0, 0.2];     % Magnitude limits (always 0 to positive)
opts.fps = fs/10;
opts.colormap = 'turbo';          % Modern alternative to jet
opts.mag_colormap = 'hot';

video_folder='G:\Shared drives\AirSeaLab\Projects\SOARS\Data\20250108\';
make_slope_video(Sx, Sy, time(500:510), strcat(video_folder,'slope_gauge_test_20260108.mp4'), opts);
%% Statistics
fprintf('\n===== SLOPE STATISTICS =====\n');
Sx_all = Sx(:);
Sy_all = Sy(:);
slope_mag_all = sqrt(Sx_all.^2 + Sy_all.^2);

fprintf('Sx: mean=%.4f, std=%.4f, range=[%.4f, %.4f]\n', ...
    mean(Sx_all), std(Sx_all), min(Sx_all), max(Sx_all));
fprintf('Sy: mean=%.4f, std=%.4f, range=[%.4f, %.4f]\n', ...
    mean(Sy_all), std(Sy_all), min(Sy_all), max(Sy_all));
fprintf('|S|: mean=%.4f, std=%.4f, max=%.4f\n', ...
    mean(slope_mag_all), std(slope_mag_all), max(slope_mag_all));

%% Time series at a point
fprintf('\n===== TIME SERIES EXAMPLE =====\n');
% Pick center point
center_y = round(height/2);
center_x = round(width/2);

figure('Name', 'Time Series at Center Point');
subplot(3,1,1);
plot(time, squeeze(Sx(center_y, center_x, :)));
xlabel('Time (s)'); ylabel('S_x');
title(sprintf('S_x at center point (%d, %d)', center_x, center_y));
grid on;

subplot(3,1,2);
plot(time, squeeze(Sy(center_y, center_x, :)));
xlabel('Time (s)'); ylabel('S_y');
title('S_y at center point');
grid on;

subplot(3,1,3);
slope_mag_center = squeeze(sqrt(Sx(center_y, center_x, :).^2 + Sy(center_y, center_x, :).^2));
plot(time, slope_mag_center);
xlabel('Time (s)'); ylabel('|S|');
title('Slope magnitude at center point');
grid on;

fprintf('\n✓ All processing complete!\n');
fprintf('Data saved in: CISG_slopes_all_frames.mat\n');
fprintf('  - Sx: [%d x %d x %d] slope array\n', height, width, totalFrames);
fprintf('  - Sy: [%d x %d x %d] slope array\n', height, width, totalFrames);
fprintf('  - time: [%d x 1] time vector (seconds)\n', totalFrames);
fprintf('  - setup: struct with experimental parameters\n');

%% Compute (ak)^2 time series
nFrames=length(1:10:end_frame);
ak2_Sx = zeros(1, nFrames);
ak2_Sy = zeros(1, nFrames);

for i = 1:nFrames
    Sx_frame = double(m.Sx(setup.roi.row_min:setup.roi.row_max, ...
        setup.roi.col_min:setup.roi.col_max, i));
    Sy_frame = double(m.Sy(setup.roi.row_min:setup.roi.row_max, ...
        setup.roi.col_min:setup.roi.col_max, i));
    ak2_Sx(i) = mean(Sx_frame(:).^2);
    ak2_Sy(i) = mean(Sy_frame(:).^2);
end

figure('Name', 'Mean Square Slope Evolution');
plot(time, ak2_Sx,'LineWidth',1.5); hold on;
plot(time, ak2_Sy,'LineWidth',1.5);
plot(time, ak2_Smag, '--k','LineWidth',1.5);
xlabel('Time (s)');
ylabel('$(ak)^2$', 'Interpreter', 'latex');
legend('S_x', 'S_y', '|S|');
grid off;
set(gca,'fontsize',12)