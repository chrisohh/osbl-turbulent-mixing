%% Generate Color Imaging Slope Gauge (CISG) Pattern
% Creates a unique RGB color at each pixel position

% Calculate pixels needed for specific print size
print_width_cm = 20*2.54;  % A3 width in cm
print_height_cm = 20*2.54; % A3 height in cm
dpi = 600;              % Print resolution

width = round(print_width_cm / 2.54 * dpi);
height = round(print_height_cm / 2.54 * dpi);

% Image dimensions (pixels) - adjust to your needs
% width = 1000;   % pixels in X direction
% height = 1000;  % pixels in Y direction

% Create coordinate matrices
[X, Y] = meshgrid(1:width, 1:height);

% Normalize coordinates to 0-255 range
Red = uint8(255 * (X - 1) / (width - 1));      % Increases left to right
Green = uint8(255 * (Y - 1) / (height - 1));   % Increases bottom to top
Blue = uint8(128 * ones(height, width));       % Constant value

% Combine into RGB image
gradient_image = cat(3, Red, Green, Blue);


% %% ADD CALIBRATION MARKERS HERE (optional)
% % White squares at corners for alignment reference
% marker_size = 20;  % size in pixels
% 
% % Top-left corner
% gradient_image(1:marker_size, 1:marker_size, 1) = 255;  % Red channel
% gradient_image(1:marker_size, 1:marker_size, 2) = 255;  % Green channel
% gradient_image(1:marker_size, 1:marker_size, 3) = 255;  % Blue channel
% 
% % Top-right corner
% gradient_image(1:marker_size, (end-marker_size+1):end, 1) = 255;
% gradient_image(1:marker_size, (end-marker_size+1):end, 2) = 255;
% gradient_image(1:marker_size, (end-marker_size+1):end, 3) = 255;
% 
% % Bottom-left corner
% gradient_image((end-marker_size+1):end, 1:marker_size, 1) = 255;
% gradient_image((end-marker_size+1):end, 1:marker_size, 2) = 255;
% gradient_image((end-marker_size+1):end, 1:marker_size, 3) = 255;
% 
% % Bottom-right corner
% gradient_image((end-marker_size+1):end, (end-marker_size+1):end, 1) = 255;
% gradient_image((end-marker_size+1):end, (end-marker_size+1):end, 2) = 255;
% gradient_image((end-marker_size+1):end, (end-marker_size+1):end, 3) = 255;

% Display the gradient
figure;
imshow(gradient_image);
title('CISG Color Gradient Pattern');
xlabel('X direction (Red increases →)');
ylabel('Y direction (Green increases ↓)');

% Save as high-quality PNG for printing
imwrite(gradient_image, 'CISG_gradient.png');

% Optional: Save as TIFF for even higher quality printing
imwrite(gradient_image, 'CISG_gradient.tif', 'Compression', 'none');

fprintf('Gradient image saved!\n');
fprintf('Image size: %d x %d pixels\n', width, height);

