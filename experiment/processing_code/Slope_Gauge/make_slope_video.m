function make_slope_video(Sx, Sy, time, x_cm, y_cm, output_filename, options)
% MAKE_SLOPE_VIDEO Create video with Sx and Sy slope fields
%
% Inputs:
%   Sx - slope field in x-direction [height x width x time]
%   Sy - slope field in y-direction [height x width x time]
%   time - time vector in seconds [1 x nframes]
%   x_cm - x-axis coordinates in cm [1 x width]
%   y_cm - y-axis coordinates in cm [1 x height]
%   output_filename - output video filename (e.g., 'slopes.mp4' or 'slopes.avi')
%   options - (optional) struct with fields:
%       .caxis_limits - [min max] color axis for Sx/Sy (default: [-0.4, 0.4])
%       .fps - frames per second (default: 30)
%       .quality - video quality 0-100 (default: 90)
%       .figure_position - [left bottom width height] (default: [100 100 1000 600])
%       .colorbar_label - label for colorbar (default: '$ak$')
%
% Example:
%   opts = struct();
%   opts.caxis_limits = [-0.4, 0.4];
%   opts.fps = 50;
%   make_slope_video(Sx, Sy, time, x_cm, y_cm, 'slopes.mp4', opts);

% Parse optional arguments
if nargin < 7
    options = struct();
end

% Set defaults
if ~isfield(options, 'caxis_limits'), options.caxis_limits = [-0.4, 0.4]; end
if ~isfield(options, 'fps'), options.fps = 30; end
if ~isfield(options, 'quality'), options.quality = 90; end
if ~isfield(options, 'figure_position'), options.figure_position = [100 100 1000 600]; end
if ~isfield(options, 'colorbar_label'), options.colorbar_label = '$ak$'; end

% Get dimensions
[height, width, nframes] = size(Sx);

fprintf('===== CREATING SLOPE VIDEO =====\n');
fprintf('Output file: %s\n', output_filename);
fprintf('Frames: %d\n', nframes);
fprintf('Frame rate: %d fps\n', options.fps);
fprintf('Duration: %.2f seconds\n', nframes/options.fps);
fprintf('Color axis: [%.3f, %.3f]\n', options.caxis_limits(1), options.caxis_limits(2));

% Create video writer
[~, ~, ext] = fileparts(output_filename);
if strcmpi(ext, '.mp4')
    v = VideoWriter(output_filename, 'MPEG-4');
    v.Quality = options.quality;
elseif strcmpi(ext, '.avi')
    v = VideoWriter(output_filename, 'Motion JPEG AVI');
    v.Quality = options.quality;
else
    error('Unsupported video format. Use .mp4 or .avi');
end
v.FrameRate = options.fps;
open(v);

% Create figure
fig = figure('Position', options.figure_position, 'Color', 'white');

% Sx plot (left panel)
subplot(1,2,1)
im1 = imagesc(x_cm, y_cm, Sx(:,:,1));
axis equal tight;
c1 = colorbar;
colormap(gca, 'redblue');
caxis(options.caxis_limits);
title(sprintf('$S_x$ at $t = %.2f$ s', time(1)), 'Interpreter', 'latex');
xlabel('$x$ (cm)', 'Interpreter', 'latex');
ylabel('$y$ (cm)', 'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'FontName', 'Times');
c1.Label.String = options.colorbar_label;
c1.Label.Interpreter = 'latex';
c1.Label.FontSize = 16;
ax1 = gca;

% Sy plot (right panel)
subplot(1,2,2)
im2 = imagesc(x_cm, y_cm, Sy(:,:,1));
axis equal tight;
c2 = colorbar;
colormap(gca, 'redblue');
caxis(options.caxis_limits);
title(sprintf('$S_y$ at $t = %.2f$ s', time(1)), 'Interpreter', 'latex');
xlabel('$x$ (cm)', 'Interpreter', 'latex');
ylabel('$y$ (cm)', 'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'FontName', 'Times');
c2.Label.String = options.colorbar_label;
c2.Label.Interpreter = 'latex';
c2.Label.FontSize = 16;
ax2 = gca;

% Write frames
fprintf('Writing frames: ');
tic;
for frame = 1:nframes
    % Update image data
    set(im1, 'CData', Sx(:,:,frame));
    set(im2, 'CData', Sy(:,:,frame));
    
    % Update titles
    set(get(ax1, 'Title'), 'String', sprintf('$S_x$ at $t = %.2f$ s', time(frame)));
    set(get(ax2, 'Title'), 'String', sprintf('$S_y$ at $t = %.2f$ s', time(frame)));
    
    % Capture frame
    drawnow;
    frame_data = getframe(fig);
    writeVideo(v, frame_data);
    
    % Progress update
    if mod(frame, 50) == 0 || frame == nframes
        elapsed = toc;
        fps_write = frame / elapsed;
        eta = (nframes - frame) / fps_write;
        fprintf('\n Frame %d/%d (%.1f%%) - %.1f fps - ETA: %.1f sec', ...
            frame, nframes, 100*frame/nframes, fps_write, eta);
    end
end
fprintf('\n');

% Close video
close(v);
close(fig);

total_time = toc;
fprintf('\n✓ Video saved: %s\n', output_filename);
fprintf('Total time: %.1f seconds\n', total_time);
fprintf('Average write speed: %.1f fps\n', nframes/total_time);

% Get file size
file_info = dir(output_filename);
fprintf('File size: %.1f MB\n', file_info.bytes/1e6);

end