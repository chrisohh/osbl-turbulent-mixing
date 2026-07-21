%% Load, Crop, Annotate, and Speed Up a .mov Video
% Requirements: MATLAB with Image Processing Toolbox and VideoReader support


%% ---- USER SETTINGS ----
inputFile  = 'G:\Shared drives\OSBL\2026-01 Pilot Experiments\Photos and videos\029_DSLR.mov';       % Path to your .mov file
outputFile = 'ML_Deepening_2x.mp4';  % Output filename

% Crop region [x, y, width, height] in pixels ([] = no crop)
cropRect = []%[100, 50, 640, 480];

fontSize  = 18;
fontColor = [1 1 1];        % White
bgColor   = [0 0 0 0.5];    % Semi-transparent black

%% ---- READ FIRST FRAME ----
vr = VideoReader(inputFile);
frame = readFrame(vr);      % grabs frame 1
t_sec = 0.0;                % first frame = t=0

%% ---- CROP ----
if ~isempty(cropRect)
    x = cropRect(1); y = cropRect(2);
    w = cropRect(3); h = cropRect(4);
    frame = frame(y:y+h-1, x:x+w-1, :);
end

%% ---- DISPLAY WITH AXES AND TIME LABEL ----
figure('Color','w');
ax = axes('Units','normalized','Position',[0.08 0.08 0.88 0.88]);
imshow(frame, 'Parent', ax);
ax.Visible  = 'on';
ax.FontSize = 12;
ax.Box      = 'on';

[fH, fW, ~] = size(frame);
xlim(ax, [0.5, fW+0.5]);
ylim(ax, [0.5, fH+0.5]);

xlabel(ax, 'x (px)','FontSize',13);
ylabel(ax, 'y (px)','FontSize',13);

timeStr = sprintf('t = %.3f s', t_sec);
text(ax, 0.02, 0.97, timeStr, ...
    'Units','normalized', ...
    'Color', fontColor, ...
    'FontSize', fontSize, ...
    'FontWeight','bold', ...
    'VerticalAlignment','top' );
