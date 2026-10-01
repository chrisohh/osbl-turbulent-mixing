%% Side-by-side video: dye deepening (Water_SURF) | angled view (Water_PIV)
% CoreView 95.  All cameras are triggered together at 50 Hz, so frame n is
% the same instant in both views.
%   Left : Water_SURF (longitudinal side view) as recorded, physical axes
%          from the CoreView 94 plate (util/lif_plate_scale.m), dye front
%          + surface overlaid (util/lif_dye_front.m).
%   Right: Water_PIV (angled view -- "PIV" is only the camera name), as
%          recorded, in pixels: the view is oblique, so one uniform scale
%          would be wrong across the frame.
% Title: time since wind start, t = camDelay + (n-1)/fs, from the hot-wire file.
%
% Only frames present in BOTH folders are written (require_both = true).
% With require_both = false, every frame in either folder is written and a
% missing view is left blank.  The dye-front time series is saved next to
% the video (<out_file>_front.mat) for the deepening analysis.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'util'));

%% ---------------- User settings ----------------
raw_root  = 'D:\HLAB_2026\Core';
run_num   = 95;
cam_L     = 'Water_Surf';
cam_R     = 'Water_PIV';

frame_start  = 1;
frame_step   = 5;          % every Nth acquired frame
frame_end    = [];         % [] = last available
require_both = true;

plate_num = 94;  plate_frm = 1;             % Water_SURF plate
sc.plate_scale = 1.0;      % <-- printed / designed size (plate was scaled down)
sc.square_mm   = 11;       % small plate: coarse 11x10 layout with 11 mm squares
                           % (big plate = same layout at the design 40 mm); overrides plate_scale
sc.offset_mm   = 25;       % plate was this far TOWARD the camera
sc.n_medium    = 1.33;
sc.lens_f_mm   = 35;       % <-- CHECK

dy.dye_frac       = 0.90;  % see lif_dye_front
dy.dye_thr        = [];
dy.smooth_mm      = 1.5;
dy.surf_window_cm = [4 10];
dy.surf_gap_cm    = 0.4;

hw_file = 'D:\HLAB_2026\hotwire\hotwire_20260923_113939.mat';
fs      = 50;

downsample = 2;
clims_L    = [350 600];    % counts
clims_R    = [];           % counts; [] = percentiles of the first right frame
pct_lims   = [1 99.5];
cmap       = 'bone';

out_file      = fullfile(raw_root, sprintf('CoreView_%d', run_num), ...
                         sprintf('CoreView_%d_surf_angled.mp4', run_num));
video_fps     = 10;
video_quality = 95;
fig_pos       = [40 40 1800 720];

%% ---------------- Frames ----------------
run_dir = fullfile(raw_root, sprintf('CoreView_%d', run_num));
fname   = @(cam, n) fullfile(run_dir, cam, sprintf('CoreView_%d_%s_%04d.raw', run_num, cam, n));
nL = frame_numbers(fullfile(run_dir, cam_L, sprintf('CoreView_%d_%s_*.raw', run_num, cam_L)));
nR = frame_numbers(fullfile(run_dir, cam_R, sprintf('CoreView_%d_%s_*.raw', run_num, cam_R)));
fprintf('%s: %d frames (%s)\n', cam_L, numel(nL), mat2str(nL([1 end])));
fprintf('%s: %d frames (%s)\n', cam_R, numel(nR), mat2str(nR([1 end])));

if require_both, avail = intersect(nL, nR); else, avail = union(nL, nR); end
if isempty(frame_end), frame_end = max(avail); end
frames = intersect(frame_start:frame_step:frame_end, avail);
if isempty(frames)
    error('No frames to write (require_both = %d). Wait for both folders to fill, or set require_both = false.', ...
          require_both);
end
nFrames = numel(frames);
fprintf('Writing %d frames (%d..%d step %d) -> %s\n', nFrames, frames(1), frames(end), frame_step, out_file);

%% ---------------- Scale, axes, time ----------------
fplate = fullfile(raw_root, sprintf('CoreView_%d', plate_num), cam_L, ...
                  sprintf('CoreView_%d_%s_%02d.raw', plate_num, cam_L, plate_frm));
[mmpp, c0] = lif_plate_scale(fplate, sc);

ny = 3072/downsample;  nx = 4096/downsample;
x_cm =  ((1:nx)*downsample - c0(1)) * mmpp / 10;
z_cm = -((1:ny)*downsample - c0(2)) * mmpp / 10;
xpx  = (1:nx)*downsample;
ypx  = (1:ny)*downsample;

H = load(hw_file);
if isfield(H, 'camStartElapsed') && ~isnan(H.camStartElapsed)
    camDelay = H.camStartElapsed - H.fanStartElapsed;
else
    camDelay = H.runConfig.DELAY_BEFORE_TRIG;
    warning('No camStartElapsed in %s -- using DELAY_BEFORE_TRIG = %g s.', hw_file, camDelay);
end
t_wind = camDelay + (frames - 1) / fs;

%% ---------------- Figure (built once, updated per frame) ----------------
fig = figure('Position', fig_pos, 'Color', 'w', 'Visible', 'off');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axL = nexttile(tl);
hL  = imagesc(axL, x_cm, z_cm, nan(ny, nx));
axis(axL, 'image'); set(axL, 'YDir', 'normal');
colormap(axL, cmap); caxis(axL, clims_L);
cb = colorbar(axL); cb.Label.String = 'counts';
hold(axL, 'on');
hF = plot(axL, x_cm, nan(1, nx), 'r-', 'LineWidth', 1.5);
hS = plot(axL, x_cm, nan(1, nx), 'c--', 'LineWidth', 1);
hold(axL, 'off');
legend(axL, {sprintf('dye front (%.0f\\%% recovery)', 100*dy.dye_frac), 'surface'}, ...
       'Interpreter', 'latex', 'Location', 'southeast');
xlabel(axL, '$x$ (cm)', 'Interpreter', 'latex');
ylabel(axL, '$z$ (cm, from board centre)', 'Interpreter', 'latex');
title(axL, 'Side view (Water\_SURF): dye deepening', 'Interpreter', 'latex');
set(axL, 'FontSize', 13, 'FontName', 'times');

axR = nexttile(tl);
hR  = imagesc(axR, xpx, ypx, nan(ny, nx));
axis(axR, 'image');                             % YDir reverse: as the camera sees it
colormap(axR, cmap);
cb = colorbar(axR); cb.Label.String = 'counts';
xlabel(axR, 'column (px)'); ylabel(axR, 'row (px)');
title(axR, 'Angled view (Water\_PIV)', 'Interpreter', 'latex');
set(axR, 'FontSize', 13, 'FontName', 'times');

hT = title(tl, '', 'Interpreter', 'latex', 'FontSize', 16);

%% ---------------- Write ----------------
v = VideoWriter(out_file, 'MPEG-4');
v.FrameRate = video_fps;
v.Quality   = video_quality;
open(v);
cleanupObj = onCleanup(@() close(v));

z_front_all = nan(nFrames, nx);
z_surf_all  = nan(nFrames, nx);
tic;
for i = 1:nFrames
    n = frames(i);

    fL = fname(cam_L, n);
    if isfile(fL)
        I = lif_load_raw(fL);
        I = I(1:downsample:end, 1:downsample:end);
        [z_front_all(i,:), z_surf_all(i,:)] = lif_dye_front(I, z_cm, mmpp*downsample, dy);
        set(hL, 'CData', I);
    else
        set(hL, 'CData', nan(ny, nx));
    end
    set(hF, 'YData', z_front_all(i,:));
    set(hS, 'YData', z_surf_all(i,:));

    fR = fname(cam_R, n);
    if isfile(fR)
        I = lif_load_raw(fR);
        I = I(1:downsample:end, 1:downsample:end);
        if isempty(clims_R), clims_R = prctile(I(:), pct_lims); caxis(axR, clims_R); end
        set(hR, 'CData', I);
    else
        set(hR, 'CData', nan(ny, nx));
    end

    hT.String = sprintf('CoreView %d, frame %d, $t = %.2f$ s since wind start', run_num, n, t_wind(i));

    rgb = print(fig, '-RGBImage', '-r0');      % works with an invisible figure
    rgb = rgb(1:2*floor(end/2), 1:2*floor(end/2), :);   % H.264 needs even size
    writeVideo(v, rgb);

    if mod(i, 25) == 0 || i == nFrames
        el = toc;
        fprintf('  %d/%d  frame %d  t = %.2f s  --  %.2f fps  ETA %.0f s\n', ...
                i, nFrames, n, t_wind(i), i/el, (nFrames-i)/(i/el));
    end
end
clear cleanupObj;   % closes the writer
fprintf('\nWrote %s  (%d frames, %.1f s)\n', out_file, nFrames, toc);

%% ---------------- Save the deepening time series ----------------
depth_all = z_surf_all - z_front_all;           % [nFrames x nx] cm below surface
[p, b] = fileparts(out_file);
ts_file = fullfile(p, [b '_front.mat']);
save(ts_file, 'frames', 't_wind', 'x_cm', 'z_front_all', 'z_surf_all', 'depth_all', ...
     'mmpp', 'c0', 'downsample', 'dy', 'sc', 'hw_file');
fprintf('Saved dye-front time series to %s\n', ts_file);

%% ---------------- helper ----------------
function n = frame_numbers(pattern)
    % Skips files still being copied (0 bytes / short), which lif_load_raw rejects
    d = dir(pattern);
    d = d([d.bytes] >= 4096*3072*2);
    n = sort(cellfun(@(s) str2double(regexp(s, '_(\d+)\.raw$', 'tokens', 'once')), {d.name}));
end
