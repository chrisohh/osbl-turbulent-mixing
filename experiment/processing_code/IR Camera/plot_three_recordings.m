%% Plot two IR camera recordings side-by-side
% Rec-000014: unstratified
% Rec-000026: stratified, 10 cm freshwater layer

addpath('C:\Program Files\FLIR Systems\sdks\file\bin\Release')

%% Recording configuration
% Each row: {filename, subplot title, mask threshold, wall temp, frame to plot, clim, temp offset}
configs = {
    'D:\HLAB_2026\IR_camera\Rec-000014.ats', 'Unstratified',        17.8, 17.3, 1166, [16.0, 16.7], 16.7;
    'D:\HLAB_2026\IR_camera\Rec-000026.ats', 'Stratified (10 cm FW)', 17.8, 17.3, 1166, [16.0, 16.7], 16.7;
};

channelWidth = 0.5;  % metres

figure('Position', [50 50 1100 520], 'Color', 'w');

for k = 1:2
    rec_file    = configs{k,1};
    rec_title   = configs{k,2};
    maskThresh  = configs{k,3};
    wallTemp    = configs{k,4};
    frameToPlot = configs{k,5};
    clim_range  = configs{k,6};
    tempOffset  = configs{k,7};

    fprintf('\n=== Processing %s ===\n', rec_title);

    %% Load recording
    v = FlirMovieReader(rec_file);
    v.unit = 'temperatureFactory';

    numFrames = v.sourceInfo.presetInfo(1).numFrames;
    height    = v.sourceInfo.imageHeight;
    width     = v.sourceInfo.imageWidth;

    %% Read all frames into memory
    data = zeros(height, width, numFrames, 'single');
    reset(v);
    for n = 1:numFrames
        data(:,:,n) = step(v);
        if mod(n, 500) == 0
            fprintf('  Loaded frame %d / %d\n', n, numFrames);
        end
    end

    frameRate = v.sourceInfo.presetInfo(1).frameRate;

    %% Background (average of first 10 frames)
    img_avg = mean(data(:,:,1:10), 3);

    %% Wall detection mask
    Frame = data(:,:, frameToPlot);
    mask  = Frame < maskThresh;

    %% Masked anomaly frame
    masked_data = (Frame - img_avg) .* mask;
    masked_data(~mask) = NaN;

    %% Detect walls via edge detection on temperature threshold
    tolerance = 0.1;
    wallMask  = abs(Frame - wallTemp) < tolerance;
    edges     = edge(wallMask, 'Canny');

    leftWall_x  = nan(height, 1);
    rightWall_x = nan(height, 1);

    for row = 1:height
        edgePixels = find(edges(row, :));
        if ~isempty(edgePixels)
            leftWall_x(row)  = edgePixels(1);
            rightWall_x(row) = edgePixels(end);
        end
    end

    leftWall_y  = (1:height)';
    rightWall_y = (1:height)';

    validLeft  = ~isnan(leftWall_x);
    validRight = ~isnan(rightWall_x);

    x_left = leftWall_x(validLeft);  y_left = leftWall_y(validLeft);
    x_right = rightWall_x(validRight); y_right = rightWall_y(validRight);

    %% Remove outliers
    med_l = median(x_left);  mad_l = median(abs(x_left  - med_l));
    med_r = median(x_right); mad_r = median(abs(x_right - med_r));

    x_left_clean  = x_left(abs(x_left  - med_l) <= 3*mad_l);
    y_left_clean  = y_left(abs(x_left  - med_l) <= 3*mad_l);
    x_right_clean = x_right(abs(x_right - med_r) <= 3*mad_r);
    y_right_clean = y_right(abs(x_right - med_r) <= 3*mad_r);

    %% Fit polynomial to walls
    polyOrder   = 2;
    p_left      = polyfit(y_left_clean,  x_left_clean,  polyOrder);
    p_right     = polyfit(y_right_clean, x_right_clean, polyOrder);

    y_eval      = 1:height;
    x_left_fit  = polyval(p_left,  y_eval);
    x_right_fit = polyval(p_right, y_eval);

    %% Build physical coordinate mapping
    [pixelY, pixelX] = meshgrid(1:width, 1:height);  % note: meshgrid swaps dims
    pixelY = pixelY';
    pixelX = pixelX';

    physicalX = zeros(height, width);
    physicalY = zeros(height, width);

    for row = 1:height
        xl = polyval(p_left,  row);
        xr = polyval(p_right, row);
        pw = xr - xl;
        ppm = pw / channelWidth;

        for col = 1:width
            physicalX(row, col) = ((col - xl) / pw - 0.5) * channelWidth;
        end

        if row == 1
            physicalY(row, :) = 0;
        else
            physicalY(row, :) = physicalY(row-1, 1) + 1/ppm;
        end
    end

    physicalY = physicalY - mean(physicalY(:));

    %% Interpolate onto regular grid
    x_regular = linspace(-channelWidth/2, channelWidth/2, width);
    y_regular = linspace(min(physicalY(:)), max(physicalY(:)), height);
    [X_regular, Y_regular] = meshgrid(x_regular, y_regular);

    frameData_corrected = griddata(physicalX, physicalY, double(masked_data), ...
        X_regular, Y_regular, 'linear');

    %% Plot
    subplot(1, 2, k);
    imagesc(y_regular*100, x_regular*100, frameData_corrected' + tempOffset);
    colormap(hot);
    c = colorbar;
    axis equal tight;
    title(sprintf('%s\n$t = %.2f$ s', rec_title, frameToPlot/frameRate), ...
        'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Along-wind (cm)', 'FontSize', 13);
    if k == 1
        ylabel('Cross-wind (cm)', 'FontSize', 13);
    end
    c.Label.String = 'Temp (°C)';
    clim(clim_range);
    set(gca, 'FontSize', 12);
end

sgtitle('IR Camera: Distortion-Corrected Surface Temperature Anomaly', ...
    'FontSize', 15, 'FontWeight', 'bold');
