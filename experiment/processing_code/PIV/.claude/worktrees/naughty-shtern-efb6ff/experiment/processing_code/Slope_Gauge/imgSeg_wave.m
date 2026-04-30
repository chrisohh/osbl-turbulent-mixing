%% 1. Load video and skip to 24 seconds
v = VideoReader('G:\My Drive\OSBL\Jan2026_test\029_DSLR.MOV');
fps = 50;
start_time = 24;

v.CurrentTime = start_time;
first_frame = readFrame(v);

%% 2. Select 4 points to define crop region
figure('Name', 'Select 4 corner points');
imshow(first_frame);
title('Click 4 corners to define crop region');
[x_pts, y_pts] = ginput(4);
plot([x_pts; x_pts(1)], [y_pts; y_pts(1)], 'r-', 'LineWidth', 2);
plot(x_pts, y_pts, 'ro', 'MarkerSize', 10);
pause(1);
close;

% Get bounding box
x_min = round(min(x_pts));
x_max = round(max(x_pts));
y_min = round(min(y_pts));
y_max = round(max(y_pts));
output_width = x_max - x_min + 1;
output_height = y_max - y_min + 1;

%% 3. Process video
v.CurrentTime = start_time;
end_time=77;%v.Duration;
num_frames = round((end_time - start_time) * fps);

cropped_frames = zeros(output_height, output_width, 3, num_frames, 'uint8');

fprintf('Reading and cropping frames...\n');
frame_idx = 0;
while hasFrame(v) && frame_idx < num_frames
    frame_idx = frame_idx + 1;
    frame = readFrame(v);
    cropped_frames(:,:,:,frame_idx) = frame(y_min:y_max, x_min:x_max, :);

    if mod(frame_idx, 100) == 0
        fprintf('Frame %d/%d\n', frame_idx, num_frames);
    end
end

cropped_frames = cropped_frames(:,:,:,1:frame_idx);

%% 4. Apply image segmentation to find reflective line
%% Process all frames with optimized gap handling
eta_matrix = zeros(frame_idx, output_width);
edge_buffer=15;

fprintf('Processing frames...\n');
for i = 1:frame_idx
    frame_sub = double(cropped_frames(:,:,:,i));

    I = rgb2gray(uint8(frame_sub));

    % Create mask to ignore edge pixels
    [rows, cols] = size(I);
    mask = true(rows, cols);
    mask(1:edge_buffer, :) = false;           % Top edge
    mask(end-edge_buffer+1:end, :) = false;   % Bottom edge
    mask(:, 1:edge_buffer) = false;           % Left edge
    mask(:, end-edge_buffer+1:end) = false;   % Right edge

    % Sobel edge detection
    BW2 = edge(I, 'Sobel',0.12);

    % Apply mask
    BW2(~mask) = false;

% Remove small connected components (dots)
min_size = 5;  % Adjust this - objects smaller than this are removed
BW2_clean = bwareaopen(BW2, min_size);

% Now extract
[rows, cols] = size(BW2_clean);
eta = zeros(1, cols);
for j = 1:cols
    column = BW2_clean(:, j);
    edge_indices = find(column);
    if ~isempty(edge_indices)
        eta(j) = mean(edge_indices);
    else
        eta(j) = NaN;
    end
end

% Fix gaps and jumps
eta = fillmissing(eta, 'linear');
eta = medfilt1(eta, 21);  % Remove jumps
eta = smoothdata(eta, 'gaussian', 50);  % Smooth

eta_matrix(i, :) = eta;

if mod(i, 100) == 0
    fprintf('Frame %d/%d\n', i, frame_idx);
end
end

% Plot

I = rgb2gray(uint8(double(cropped_frames(:,:,:,i))));
figure
imshow(I,[]);hold on
plot(eta_matrix(i,:),'r','LineWidth',2)
set(gca,'YDir','reverse')

%%
%% Create video with detected surface line overlay
output_filename = 'surface_detection_video.mp4';
video_fps = 50;  % Output video frame rate (can be different from original)

% Setup video writer
v_out = VideoWriter(output_filename, 'MPEG-4');
v_out.FrameRate = video_fps;
v_out.Quality = 95;
open(v_out);

% Create figure for rendering
fig = figure('Position', [100 100 output_width output_height], ...
    'Color', 'k');
ax = axes('Position', [0 0 1 1]);

fprintf('Creating video...\n');


for i = 1:frame_idx
    % Get grayscale frame
    I = rgb2gray(uint8(double(cropped_frames(:,:,:,i))));

    % Display frame
    imshow(I, []);
    hold on;

    % Overlay detected surface line
    plot(eta_matrix(i,:), 'r', 'LineWidth', 2);

    % Optional: Add frame number
    text(10, 20, sprintf('Frame: %d', i), ...
        'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold');

    set(gca, 'YDir', 'reverse');
    hold off;

    % Capture frame and write to video
    frame_out = getframe(ax);
    writeVideo(v_out, frame_out);

    if mod(i, 100) == 0
        fprintf('Rendered %d/%d frames\n', i, frame_idx);
    end
end

close(v_out);
close(fig);

fprintf('Video saved as: %s\n', output_filename);