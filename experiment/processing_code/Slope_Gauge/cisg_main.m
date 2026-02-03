%% CISG Main Script - Complete Workflow for CoreView Data
clear all;

%% Setup parameters
ref_folderName = "D:\CoreView_2\Flare 12M125 CCL C1112A00\";
obs_folderName = ref_folderName;
fs = 50;  % Hz
dt = 1/fs;  % 0.02 seconds between frames

% Physical setup (ADJUST THESE TO YOUR ACTUAL MEASUREMENTS)
setup = struct();
setup.camera_height = 60;            % cm above water surface (MEASURE THIS)
setup.water_depth = 50;              % cm from water surface to pattern (MEASURE THIS)
setup.pattern_width_cm = 20*2.54;    % 50.8 cm (your printed pattern)
setup.pattern_height_cm = 20*2.54;   % 50.8 cm
setup.n_water = 1.333;
setup.n_air = 1.000;
setup.frame_rate = fs;
setup.dt = dt;

fprintf('===== CISG PROCESSING SETUP =====\n');
fprintf('Data folder: %s\n', ref_folderName);
fprintf('Frame rate: %d Hz (dt = %.4f s)\n', fs, dt);
fprintf('Camera height: %.1f cm\n', setup.camera_height);
fprintf('Water depth: %.1f cm\n', setup.water_depth);
fprintf('Pattern size: %.1f x %.1f cm\n\n', setup.pattern_width_cm, setup.pattern_height_cm);

%% Load and display reference image
ref_fileName = "CoreView_2_Flare 12M125 CCL C1112A00_001.raw";
fprintf('Loading reference image: %s\n', ref_fileName);
[ref_image, meta_ref] = cisg_load_coreview(strcat(ref_folderName, ref_fileName));

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

%% Process all frames
totalFrames = 753;
fprintf('===== PROCESSING %d FRAMES =====\n', totalFrames);

% Time vector
time = (0:totalFrames-1) * dt;

% Process each frame
tic;
count=0;
for n = 500:510;%1:totalFrames
    count=1+count;
    % Load observed image
    obs_fileName = sprintf("CoreView_2_Flare 12M125 CCL C1112A00_%03d.raw", n);
    [obs_image, meta_obs] = cisg_load_coreview(strcat(obs_folderName, obs_fileName));
    
    % Calculate slopes
    [Sx(:,:,count), Sy(:,:,count)] = cisg_calculate_slopes(ref_image, obs_image, setup);
    
    % Progress update every 50 frames
    if mod(n, 50) == 0
        elapsed = toc;
        fps_processing = n / elapsed;
        est_remaining = (totalFrames - n) / fps_processing;
        fprintf('  Frame %d/%d (%.1f%%) - %.1f fps - ETA: %.1f sec\n', ...
            n, totalFrames, 100*n/totalFrames, fps_processing, est_remaining);
    end
end
total_time = toc;

fprintf('\n✓ Processing complete!\n');
fprintf('  Total time: %.1f seconds\n', total_time);
fprintf('  Average speed: %.1f fps\n', totalFrames/total_time);

%% Save results
fprintf('\n===== SAVING RESULTS =====\n');
save('CISG_slopes_all_frames.mat', 'Sx', 'Sy', 'time', 'setup', '-v7.3');
fprintf('Saved: CISG_slopes_all_frames.mat\n');

% Save metadata
fid = fopen('CISG_processing_info.txt', 'w');
fprintf(fid, 'CISG Processing Summary\n');
fprintf(fid, '======================\n\n');
fprintf(fid, 'Data folder: %s\n', ref_folderName);
fprintf(fid, 'Reference file: %s\n', ref_fileName);
fprintf(fid, 'Total frames: %d\n', totalFrames);
fprintf(fid, 'Frame rate: %d Hz\n', fs);
fprintf(fid, 'Duration: %.2f seconds\n', totalFrames * dt);
fprintf(fid, 'Image size: %d x %d\n', height, width);
fprintf(fid, 'Camera height: %.1f cm\n', setup.camera_height);
fprintf(fid, 'Water depth: %.1f cm\n', setup.water_depth);
fprintf(fid, 'Pattern size: %.1f x %.1f cm\n', setup.pattern_width_cm, setup.pattern_height_cm);
fprintf(fid, 'Processing time: %.1f seconds\n', total_time);
fclose(fid);
fprintf('Saved: CISG_processing_info.txt\n');

%% Quick visualization of results
fprintf('\n===== VISUALIZING SAMPLE FRAMES =====\n');

figure('Name', 'Spatial Maps', 'Position', [100 100 1600 1000]);
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

frames_to_show = [1, round(size(Sx,3)/2), size(Sx,3)];

for i = 1:length(frames_to_show)
    frame_num = frames_to_show(i);
    
    nexttile;
    imagesc(Sx(:,:,i));
    axis equal tight; colorbar;
    title(sprintf('S_x: Frame %d (t=%.2fs)', frame_num, time(frame_num)));
    colormap(gca, 'turbo');
    caxis([-0.2,0.2]);
end

for i = 1:length(frames_to_show)
    frame_num = frames_to_show(i);
    
    nexttile;
    imagesc(Sy(:,:,i));
    axis equal tight; colorbar;
    title(sprintf('S_y: Frame %d', frame_num));
    colormap(gca, 'turbo');
    caxis([-0.2,0.2]);
end


for i = 1:length(frames_to_show)
    frame_num = frames_to_show(i);
    Smag = sqrt(Sx(:,:,i).^2 + Sy(:,:,i).^2);
    nexttile;
    imagesc(Smag);
    axis equal tight; colorbar;
    title(sprintf('|S|: Frame %d', frame_num));
    colormap(gca, 'hot');
    caxis([0,0.2]);
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