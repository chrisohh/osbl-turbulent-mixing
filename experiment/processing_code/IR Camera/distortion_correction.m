
%% Directly read from ats
addpath('C:\Program Files\FLIR Systems\sdks\file\bin\Release')
% Create a FlirMovieReader object
v = FlirMovieReader('\\Airseaserver28\D\HLAB_2026\IR_camera\Rec-000014.ats');
v.unit='temperatureFactory';
%% Extract key metadata
% Camera and acquisition info
cameraInfo.model = v.sourceInfo.cameraModel;
cameraInfo.serial = v.sourceInfo.cameraSerial;
cameraInfo.lens = v.sourceInfo.lens;
cameraInfo.width = v.sourceInfo.imageWidth;
cameraInfo.height = v.sourceInfo.imageHeight;
cameraInfo.frameRate = v.sourceInfo.presetInfo(1).frameRate;
cameraInfo.numFrames = v.sourceInfo.presetInfo(1).numFrames;
cameraInfo.integrationTime = v.sourceInfo.presetInfo(1).intTime;

%% Preallocate 3D matrix for all frames
numFrames = cameraInfo.numFrames;
height = cameraInfo.height;
width = cameraInfo.width;

data = zeros(height, width, numFrames, 'single');

%% Read one frame
sampleFrame=readFrame(v, 8514);
figure;
imagesc(sampleFrame)
colorbar;
colormap(hot);
title(sprintf('t= %.2f s', frameToPlot/50));
%% Define frame range to load
startFrameIdx = 1;
endFrameIdx = numFrames;  % or specify: endFrameIdx = 1000;

numFramesToLoad = endFrameIdx - startFrameIdx + 1;
%% Skip to start frame
reset(v);
fprintf('Skipping to frame %d...\n', startFrameIdx);
for i = 1:(startFrameIdx-1)
    step(v);
end

%% Read frames
fprintf('Reading frames...\n');
frameIdx = 1;
actualFrameNum = startFrameIdx;

while ~isDone(v) && actualFrameNum <= endFrameIdx
    data(:,:,frameIdx) = step(v);
    
    if mod(frameIdx, 100) == 0
        fprintf('  Loaded frame %d/%d (source frame %d)\n', ...
                frameIdx, numFramesToLoad, actualFrameNum);
    end
    
    frameIdx = frameIdx + 1;
    actualFrameNum = actualFrameNum + 1;
end

%% Plot a sample frame
frameToPlot = 1166; % or round(numFrames/2) for middle frame

% figure('Position', [100 100 1200 500]);

imagesc(data(:,:,frameToPlot));
axis image;
colorbar;
colormap(hot);
% caxis([radiometry.minTemp radiometry.maxTemp]); % or use percentiles
title(sprintf('t= %.2f s', frameToPlot/50));
xlabel('X (pixels)');
ylabel('Y (pixels)');

%%
% % Initialize parameters
% file_numbers = 762:772;
% num_files = length(file_numbers);
% file_path='\\Airseaserver28\D\HLAB_2026\IR_camera\Rec-000014_background\';
% % Load first file to get dimensions
% first_file = sprintf('Rec-000014_background_%d.mat', file_numbers(1));
% data = load(strcat(file_path,first_file));
% img_sum = double(data.Frame);  % Convert to double for accuracy
% 
% % Loop through remaining files and accumulate
% for i = 2:num_files
%     filename = sprintf('Rec-000014_background_%d.mat', file_numbers(i));
%     data = load(strcat(file_path,filename));
%     img_sum = img_sum + double(data.Frame);
% end
% 
% % Compute time average
% img_avg = img_sum / num_files;
% 
% % Display result
% figure;
% imagesc(img_avg);
% colorbar;
% title('Time-Averaged Image');
% axis equal tight;

%%
img_sum=zeros(size(data(:,:,1)));
for i = 1:10
    img_sum = img_sum + data(:,:,i);
end

% Compute time average
img_avg = img_sum / 10;

%%
f1=data(:,:,1166);
mask = f1 < 17.8;

%%
Frame=data(:,:,1166);
masked_data=(Frame-img_avg).*mask;
masked_data(~mask) = NaN;  % Set false regions to NaN
%%

figure;
imagesc(masked_data);
colormap('hot')
colorbar
xlabel('')

%%
% Threshold to find walls (adjust based on your image)
wallTemp = 17.3;  % or whatever temp the walls are at
tolerance = 0.1;  % degrees

wallMask = abs(Frame - wallTemp) < tolerance;

figure;
imagesc(wallMask);
colormap(gray);
title('Wall detection mask');

% Find edges
edges = edge(wallMask, 'Canny');

% Extract left and right edges
% This is a simplified version - you may need to refine
figure;
imagesc(edges);
title('Detected edges - verify walls are captured');

%% Step 2: Identify left and right walls
% Assume walls are the leftmost and rightmost continuous edges
% For each row, find the leftmost and rightmost edge pixels

leftWall_x = nan(512, 1);
leftWall_y = (1:512)';
rightWall_x = nan(512, 1);
rightWall_y = (1:512)';

for row = 1:512
    edgePixels = find(edges(row, :));
    
    if ~isempty(edgePixels)
        % Left wall is the leftmost edge in this row
        leftWall_x(row) = edgePixels(1);
        % Right wall is the rightmost edge in this row
        rightWall_x(row) = edgePixels(end);
    end
end

% Remove NaN values (rows where no edges were detected)
validLeft = ~isnan(leftWall_x);
validRight = ~isnan(rightWall_x);

x_left = leftWall_x(validLeft);
y_left = leftWall_y(validLeft);
x_right = rightWall_x(validRight);
y_right = rightWall_y(validRight);

% Plot detected wall points
subplot(2,3,3);
imagesc(Frame);
colormap(hot);
colorbar;
hold on;
plot(x_left, y_left, 'g.', 'MarkerSize', 3);
plot(x_right, y_right, 'b.', 'MarkerSize', 3);
axis equal tight;
title('Detected wall points');
legend('Left wall', 'Right wall');

%% Step 3: Filter outliers using robust fitting
% Remove obvious outliers before fitting

% For left wall - remove points too far from median
median_left = median(x_left);
mad_left = median(abs(x_left - median_left)); % Median absolute deviation
outliers_left = abs(x_left - median_left) > 3 * mad_left;
x_left_clean = x_left(~outliers_left);
y_left_clean = y_left(~outliers_left);

% For right wall
median_right = median(x_right);
mad_right = median(abs(x_right - median_right)); 
outliers_right = abs(x_right - median_right) > 3 * mad_right;
x_right_clean = x_right(~outliers_right);
y_right_clean = y_right(~outliers_right);

fprintf('Left wall: %d points, %d outliers removed\n', length(x_left), sum(outliers_left));
fprintf('Right wall: %d points, %d outliers removed\n', length(x_right), sum(outliers_right));

%% Step 4: Fit polynomial curves to walls
% Higher polynomial order captures barrel/pincushion distortion
polyOrder = 2;  % Start with 5th order, adjust if needed

% Fit left wall: x as function of y
p_left = polyfit(y_left_clean, x_left_clean, polyOrder);

% Fit right wall: x as function of y
p_right = polyfit(y_right_clean, x_right_clean, polyOrder);

% Evaluate fitted curves
y_eval = 1:512;
x_left_fit = polyval(p_left, y_eval);
x_right_fit = polyval(p_right, y_eval);

% Plot fitted curves
subplot(2,3,4);
imagesc(Frame);
colormap(hot);
colorbar;
hold on;
plot(x_left_clean, y_left_clean, 'g.', 'MarkerSize', 3);
plot(x_right_clean, y_right_clean, 'b.', 'MarkerSize', 3);
plot(x_left_fit, y_eval, 'g-', 'LineWidth', 3);
plot(x_right_fit, y_eval, 'b-', 'LineWidth', 3);
axis equal tight;
title(sprintf('Fitted walls (order %d)', polyOrder));
legend('Left points', 'Right points', 'Left fit', 'Right fit');


%% Step 5: Analyze distortion
% Plot wall separation vs. y position to see distortion pattern
wall_separation = x_right_fit - x_left_fit;
mean_separation = mean(wall_separation);

%% Step 6: Create undistorted coordinate mapping
% % Known channel width
% channelWidth = 0.5;  % meters
% 
% [pixelY, pixelX] = meshgrid(1:512, 1:640);
% pixelY = pixelY';
% pixelX = pixelX';
% 
% % Initialize physical coordinate arrays
% physicalX = zeros(size(pixelX));
% physicalY = zeros(size(pixelY));
% 
% % For each row (along-wind position)
% for row = 1:512
%     % Get wall positions at this row
%     x_left_row = polyval(p_left, row);
%     x_right_row = polyval(p_right, row);
% 
%     % Current pixel width between walls
%     pixel_width = x_right_row - x_left_row;
% 
%     % Map pixel positions to physical cross-wind positions
%     % Left wall = -0.25 m, Right wall = +0.25 m
%     for col = 1:640
%         % Normalize position between walls (0 to 1)
%         normalized_pos = (col - x_left_row) / pixel_width;
% 
%         % Map to physical cross-wind position
%         physicalX(row, col) = (normalized_pos - 0.5) * channelWidth;
%     end
% end
% 
% % Scale Y direction
% % Use average pixel spacing across channel
% pixelsPerMeter_x = mean_separation / channelWidth;
% 
% % Assume square pixels in physical space
% pixelsPerMeter_y = pixelsPerMeter_x;
% physicalY = (pixelY - 256) / pixelsPerMeter_y;  % Center at middle row
% 
% fprintf('\nSpatial calibration:\n');
% fprintf('Average pixels across channel: %.2f pixels\n', mean_separation);
% fprintf('X resolution: %.3f mm/pixel\n', 1000 * channelWidth / mean_separation);
% fprintf('Y resolution: %.3f mm/pixel (assuming square pixels)\n', 1000 / pixelsPerMeter_y);
%% Step 6: Create undistorted coordinate mapping using wall geometry
% Known channel width
channelWidth = 0.5;  % meters

[pixelY, pixelX] = meshgrid(1:512, 1:640);
pixelY = pixelY';
pixelX = pixelX';

% Initialize physical coordinate arrays
physicalX = zeros(size(pixelX));
physicalY = zeros(size(pixelY));

% Calculate physical Y position from wall separation
% Key insight: wall separation in pixels tells us about perspective scaling
for row = 1:512
    % Get wall positions at this row
    x_left_row = polyval(p_left, row);
    x_right_row = polyval(p_right, row);

    % Current pixel width between walls
    pixel_width = x_right_row - x_left_row;

    % Local scale factor (pixels per meter) in X direction at this row
    pixelsPerMeter_local = pixel_width / channelWidth;

    % Map pixel positions to physical cross-wind positions
    for col = 1:640
        % Normalize position between walls (0 to 1)
        normalized_pos = (col - x_left_row) / pixel_width;

        % Map to physical cross-wind position
        physicalX(row, col) = (normalized_pos - 0.5) * channelWidth;
    end

    % Calculate physical Y position for this row
    % Integrate the local pixel spacing
    if row == 1
        physicalY(row, :) = 0;  % Reference row
    else
        % Distance between rows in physical space
        % dy_physical = dy_pixel / pixelsPerMeter_local
        dy_pixel = 1;  % One pixel step
        dy_physical = dy_pixel / pixelsPerMeter_local;
        physicalY(row, :) = physicalY(row-1, 1) + dy_physical;
    end
end

% Center Y coordinate
physicalY = physicalY - mean(physicalY(:));


%% Step 7: Create regular grid and interpolate
% Create regular physical grid
x_regular = linspace(-channelWidth/2, channelWidth/2, 640);  % -0.25 to +0.25 m
y_range = max(physicalY(:)) - min(physicalY(:));
y_regular = linspace(min(physicalY(:)), max(physicalY(:)), 512);

[X_regular, Y_regular] = meshgrid(x_regular, y_regular);

fprintf('\nPhysical domain after correction:\n');
fprintf('X: %.3f m (%.1f cm) - cross-wind\n', channelWidth, channelWidth*100);
fprintf('Y: %.3f m (%.1f cm) - along-wind\n', y_range, y_range*100);
fprintf('Aspect ratio: %.3f\n', y_range / channelWidth);

% Interpolate onto regular grid
frameData_corrected = griddata(physicalX, physicalY, double(masked_data), ...
    X_regular, Y_regular, 'linear');


%% Step 8: Visualize correction
figure('Position', [100 100 1400 600]);

subplot(1,2,1);
imagesc(masked_data);
colormap(hot);
colorbar;
axis equal tight;
hold on;
plot(x_left_fit, y_eval, 'b-', 'LineWidth', 2);
plot(x_right_fit, y_eval, 'b-', 'LineWidth', 2);
title('Original');
xlabel('Pixel X');
ylabel('Pixel Y');
legend('Wall outline');
% clim([-0.2,0.2])
subplot(1,2,2);
imagesc(x_regular*100, y_regular*100, frameData_corrected);
colormap(hot);
colorbar;
axis equal tight;
title('Distortion corrected');
xlabel('Cross-wind (cm)');
ylabel('Along-wind (cm)');
% clim([-0.2,0.2])

%% Final image
figure;
imagesc(y_regular*100, x_regular*100, frameData_corrected'+16.7);
colormap(hot);
c=colorbar;
axis equal tight;
title(sprintf('$t=%.2f$ s',1166/cameraInfo.frameRate),'interpreter','latex');
xlabel('Along-wind (cm)')
ylabel('Cross-wind (cm)');
c.Label.String='Temp (^oC)';
set(gca,'fontsize',16)
clim([16,16.7])
%% Do it for the rest of the frames
outputVideo = VideoWriter('G:\Shared drives\OSBL\2026-01 Pilot Experiments\Processed_videos\unstrat_IR_Rec_000014.mp4', 'MPEG-4');
outputVideo.FrameRate = 10;
outputVideo.Quality = 100;
open(outputVideo);

fig = figure('Visible', 'off', 'Position', [100 100 800 600]);
colormap(hot);

fprintf('Processing %d frames...\n', size(data,3));

for n = 1000:5:4000%size(data,3)
    if mod(n, 100) == 0
        fprintf('Frame %d / %d\n', n, size(data,3));
    end
    
    masked_data = (data(:,:,n) - img_avg) .* mask;
    masked_data(~mask) = NaN;
    frameData_corrected = griddata(physicalX, physicalY, double(masked_data), ...
        X_regular, Y_regular, 'linear');
    
    imagesc(y_regular*100, x_regular*100, frameData_corrected'+16.7);
    colormap(hot);
    c = colorbar;
    axis equal tight;
    title(sprintf('$t=%.2f$ s', n/cameraInfo.frameRate), 'interpreter', 'latex');
    xlabel('Along-wind (cm)');
    ylabel('Cross-wind (cm)');
    c.Label.String = 'Temp (^oC)'
    set(gca, 'fontsize', 16);
    clim([16, 16.7]);
    
    frame = getframe(fig);
    writeVideo(outputVideo, frame);
end

close(outputVideo);
close(fig);
fprintf('Video saved as corrected_frames.mp4\n');