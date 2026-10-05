%% Water_PIV video: angled view, time in the title
% One panel, one video frame per camera frame: the Water_PIV angled view as
% recorded, in pixels ("PIV" is only the camera name).  The view is oblique,
% so one uniform cm scale would be wrong across the frame -- no physical
% axes, no surface, no z99, and no hot-wire panel.
% Frame/time bookkeeping is the same as theta_hw_video.m: the title carries
% t = camDelay + (n-1)/fs, s since wind start, with camDelay from the
% hot-wire file (or its cached sync).  Set hw_file = [] to time from the
% camera start instead, t = (n-1)/fs.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'C:\Users\Administrator\Downloads\';
cam       = 'Water_PIV';
run_num   = 95;

frame_start = 1432;
frame_step  = 10;          % every Nth acquired frame
frame_end   = 1900;        % [] = last available

hw_file = 'C:\Users\Administrator\Downloads\hotwire_20260923_113939.mat';
                           % camDelay only (no velocity panel); [] = time from camera start
fs      = 50;

downsample = 2;

cmap     = 'bone';
clims    = [];             % counts; [] = percentiles of the first written frame
pct_lims = [1 99.5];
% The image keeps its true aspect (axis image), so the figure wants height
% ~ width*3072/4096 plus room for the title and labels.
fig_pos  = [40 40 1200 1000];
visible  = 'on';

out_file      = fullfile(raw_root, sprintf('CoreView_%d', run_num), ...
                         sprintf('CoreView_%d_angled.mp4', run_num));
video_fps     = 10;
video_quality = 95;

%% ---------------- Frames ----------------
cam_dir = fullfile(raw_root, sprintf('CoreView_%d', run_num), cam);
fname   = @(n) fullfile(cam_dir, sprintf('CoreView_%d_%s_%04d.raw', run_num, cam, n));
d = dir(fullfile(cam_dir, sprintf('CoreView_%d_%s_*.raw', run_num, cam)));
d = d([d.bytes] >= 4096*3072*2);                 % skip files still being copied
avail = sort(cellfun(@(s) str2double(regexp(s, '_(\d+)\.raw$', 'tokens', 'once')), {d.name}));
if isempty(avail), error('No complete frames in %s', cam_dir); end
if isempty(frame_end), frame_end = max(avail); end
frames  = intersect(frame_start:frame_step:frame_end, avail);
if isempty(frames)
    error('No frames in %d:%d:%d exist in %s (available %d..%d).', ...
          frame_start, frame_step, frame_end, cam_dir, avail(1), avail(end));
end
nFrames = numel(frames);
fprintf('%s: %d complete frames (%d..%d); writing %d (%d..%d step %d) -> %s\n', ...
        cam, numel(avail), avail(1), avail(end), nFrames, ...
        frames(1), frames(end), frame_step, out_file);

%% ---------------- Axes (pixels) and time ----------------
ny = 3072/downsample;  nx = 4096/downsample;
xpx = (1:nx)*downsample;
ypx = (1:ny)*downsample;

cal_dir = fullfile(fileparts(mfilename('fullpath')), 'calibration');   % shared caches
if isempty(hw_file)
    camDelay = 0;  t_label = 'since camera start';
else
    [~, hw_name] = fileparts(hw_file);
    camDelay = lif_hw_cached(hw_file, 0.02, fullfile(cal_dir, ['hwsync_' hw_name '.mat']));
    t_label  = 'since wind start';
end
t_wind = camDelay + (frames - 1) / fs;

%% ---------------- Figure (built once, updated per frame) ----------------
fig = figure('Position', fig_pos, 'Color', 'w', 'Visible', visible);
axI = axes(fig);
hI  = imagesc(axI, xpx, ypx, nan(ny, nx));
axis(axI, 'image');                              % YDir reverse: as the camera sees it
colormap(axI, cmap);
cb = colorbar(axI); cb.Label.String = 'counts';
xlabel(axI, 'column (px)'); ylabel(axI, 'row (px)');
set(axI, 'FontSize', 13, 'FontName', 'times');
hT = title(axI, '', 'Interpreter', 'latex', 'FontSize', 15);

%% ---------------- Writer ----------------
v = VideoWriter(out_file, 'MPEG-4');
v.FrameRate = video_fps;
v.Quality   = video_quality;
open(v);
cleanupObj = onCleanup(@() close(v));

%% ---------------- Frames ----------------
tic;
for i = 1:nFrames
    n = frames(i);
    I = lif_load_raw(fname(n));
    I = I(1:downsample:end, 1:downsample:end);

    if isempty(clims), clims = prctile(I(:), pct_lims); caxis(axI, clims); end
    set(hI, 'CData', I);
    hT.String = sprintf('CoreView %d, Water\\_PIV, frame %d, $t = %.2f$ s %s', ...
                        run_num, n, t_wind(i), t_label);

    rgb = print(fig, '-RGBImage', '-r0');                   % works with an invisible figure
    if i == 1
        vsz = 2*floor([size(rgb,1) size(rgb,2)]/2);         % H.264 needs even size
    end
    rgb = lif_fit_frame(rgb, vsz);                          % print size can drift by a pixel
    writeVideo(v, rgb);

    if mod(i, 25) == 0 || i == nFrames
        el = toc;
        fprintf('  %d/%d  frame %d  t = %.2f s  --  %.2f fps  ETA %.0f s\n', ...
                i, nFrames, n, t_wind(i), i/el, (nFrames-i)/(i/el));
    end
end
clear cleanupObj;   % closes the writer
fprintf('\nWrote %s  (%d frames, %.1f s)\n', out_file, nFrames, toc);
