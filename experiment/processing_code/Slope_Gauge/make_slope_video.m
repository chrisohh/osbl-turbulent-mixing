function make_slope_video(Sx, Sy, time, output_filename, options)
% MAKE_SLOPE_VIDEO_3PANEL Create video with Sx, Sy, and magnitude
%
% Inputs:
%   Sx - slope field in x-direction [height x width x time]
%   Sy - slope field in y-direction [height x width x time]
%   time - time vector in seconds [1 x nframes]
%   output_filename - output video filename (e.g., 'slopes.mp4' or 'slopes.avi')
%   options - (optional) struct with fields:
%     .caxis_limits - [min max] color axis for Sx/Sy (default: auto)
%     .mag_caxis_limits - [min max] color axis for magnitude (default: [0 max])
%     .fps - frames per second (default: 30)
%     .quality - video quality 0-100 (default: 90)
%     .colormap - colormap for Sx/Sy (default: 'jet')
%     .mag_colormap - colormap for magnitude (default: 'hot')
%     .figure_position - [left bottom width height] (default: [50 50 2000 900])
%     .show_colorbar - true/false (default: true)
%     .title_prefix - prefix for title (default: 'Slope Fields')
%
% Example:
%   opts = struct();
%   opts.caxis_limits = [-0.2, 0.2];
%   opts.mag_caxis_limits = [0, 0.3];
%   opts.fps = 50;
%   make_slope_video_3panel(Sx, Sy, time, 'slopes_3panel.mp4', opts);

    % Parse optional arguments
    if nargin < 5
        options = struct();
    end
    
    % Set defaults
    if ~isfield(options, 'caxis_limits') || isempty(options.caxis_limits)
        % Auto-scale based on 99th percentile
        Sx_lim = quantile(abs(Sx(:)), 0.99);
        Sy_lim = quantile(abs(Sy(:)), 0.99);
        max_lim = max(Sx_lim, Sy_lim);
        options.caxis_limits = [-max_lim, max_lim];
    end
    if ~isfield(options, 'mag_caxis_limits') || isempty(options.mag_caxis_limits)
        % Auto-scale magnitude from 0 to 99th percentile
        slope_mag = sqrt(Sx.^2 + Sy.^2);
        mag_max = quantile(slope_mag(:), 0.99);
        options.mag_caxis_limits = [0, mag_max];
    end
    if ~isfield(options, 'fps'), options.fps = 30; end
    if ~isfield(options, 'quality'), options.quality = 90; end
    if ~isfield(options, 'colormap'), options.colormap = 'jet'; end
    if ~isfield(options, 'mag_colormap'), options.mag_colormap = 'hot'; end
    if ~isfield(options, 'figure_position'), options.figure_position = [50 50 2000 900]; end
    if ~isfield(options, 'show_colorbar'), options.show_colorbar = true; end
    if ~isfield(options, 'title_prefix'), options.title_prefix = 'Slope Fields'; end
    
    % Get dimensions
    [height, width, nframes] = size(Sx);
    
    fprintf('===== CREATING 3-PANEL SLOPE VIDEO =====\n');
    fprintf('Output file: %s\n', output_filename);
    fprintf('Frames: %d\n', nframes);
    fprintf('Frame rate: %d fps\n', options.fps);
    fprintf('Duration: %.2f seconds\n', nframes/options.fps);
    fprintf('Sx/Sy color axis: [%.3f, %.3f]\n', options.caxis_limits(1), options.caxis_limits(2));
    fprintf('Magnitude color axis: [%.3f, %.3f]\n', options.mag_caxis_limits(1), options.mag_caxis_limits(2));
    
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
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Compute magnitude for first frame
    slope_mag_1 = sqrt(Sx(:,:,1).^2 + Sy(:,:,1).^2);
    
    % Tile 1: Sx
    ax1 = nexttile;
    im1 = imagesc(Sx(:,:,1));
    axis equal tight;
    if options.show_colorbar
        colorbar;
    end
    title(sprintf('S_x: t = %.3f s', time(1)));
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    colormap(ax1, options.colormap);
    caxis(options.caxis_limits);
    
    % Tile 2: Sy
    ax2 = nexttile;
    im2 = imagesc(Sy(:,:,1));
    axis equal tight;
    if options.show_colorbar
        colorbar;
    end
    title(sprintf('S_y: t = %.3f s', time(1)));
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    colormap(ax2, options.colormap);
    caxis(options.caxis_limits);
    
    % Tile 3: Magnitude
    ax3 = nexttile;
    im3 = imagesc(slope_mag_1);
    axis equal tight;
    if options.show_colorbar
        colorbar;
    end
    title(sprintf('|S|: t = %.3f s', time(1)));
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
    colormap(ax3, options.mag_colormap);
    caxis(options.mag_caxis_limits);
    
    % Overall title
    title(t, sprintf('%s - Frame 1/%d', options.title_prefix, nframes), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Write frames
    fprintf('Writing frames: ');
    tic;
    for frame = 1:nframes
        % Compute magnitude
        slope_mag = sqrt(Sx(:,:,frame).^2 + Sy(:,:,frame).^2);
        
        % Update image data
        set(im1, 'CData', Sx(:,:,frame));
        set(im2, 'CData', Sy(:,:,frame));
        set(im3, 'CData', slope_mag);
        
        % Update titles
        title(ax1, sprintf('S_x: t = %.3f s', time(frame)),'Color', 'black');
        title(ax2, sprintf('S_y: t = %.3f s', time(frame)),'Color', 'black');
        title(ax3, sprintf('|S|: t = %.3f s', time(frame)),'Color', 'black');
        title(t, sprintf('%s - Frame %d/%d', options.title_prefix, frame, nframes), ...
            'FontSize', 14, 'FontWeight', 'bold','Color', 'black');
        
        set(ax1, 'XColor', 'black', 'YColor', 'black', 'FontSize', 12);
        set(ax2, 'XColor', 'black', 'YColor', 'black', 'FontSize', 12);
        set(ax3, 'XColor', 'black', 'YColor', 'black', 'FontSize', 12);

        % Capture frame
        drawnow;
        frame_data = getframe(fig);
        writeVideo(v, frame_data);
        
        % Progress update
        if mod(frame, 50) == 0 || frame == nframes
            elapsed = toc;
            fps_write = frame / elapsed;
            eta = (nframes - frame) / fps_write;
            fprintf('\n  Frame %d/%d (%.1f%%) - %.1f fps - ETA: %.1f sec', ...
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